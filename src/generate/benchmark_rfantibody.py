#!/usr/bin/env python3
"""Prepare and run the 21-target RFantibody smoke benchmark in WSL.

Two backbones per target/arm, one sequence each, one RF2 recycle.
Resume with the same command. Failed attempts are retained and retried only
with --retry-failed. No selection thresholds are fitted to this test set.
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time

from pilot_rfantibody import inspect_pdb

ROOT = Path(__file__).resolve().parents[2]
PROC = ROOT / 'data/processed'


def stamp():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def table(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))


def prepare():
    """Validate the complete mapping before writing antigen-only inputs."""
    meta = table(PROC / 'targets.csv')
    qc = {r['PDB_ID']: r for r in table(PROC / 'structural_qc_final.csv')}
    key = lambda r: (r['PDB_ID'], r['original_chain'], r['original_res_id'],
                     r['original_insertion_code'].strip())
    mapping = {}
    for r in table(PROC / 'antigen_residue_mapping_v2.csv'):
        if key(r) in mapping or r['output_chain'] != 'T':
            raise ValueError('Duplicate mapping or unexpected target chain')
        mapping[key(r)] = int(r['output_res_id'])
    epitopes = {}
    for r in table(PROC / 'epitope_residues_v2.csv'):
        if key(r) not in mapping:
            raise ValueError(f'Unmapped epitope: {key(r)}')
        epitopes.setdefault(r['PDB_ID'], []).append(
            (float(r['min_distance_to_nanobody_A']), mapping[key(r)]))
    if len(meta) != 21 or set(epitopes) != {r['PDB_ID'] for r in meta}:
        raise ValueError('Expected exactly 21 targets with epitopes')
    prepared = []
    for r in meta:
        pid = r['PDB_ID']
        source = PROC / 'prepped_v2' / f'{pid}_H_T.pdb'
        atoms, modifications = [], []
        for s in source.read_text().splitlines():
            if not s.startswith(('ATOM', 'HETATM')) or s[21] != 'T':
                continue
            if s[17:20] == 'MSE':
                # Model parsers expect standard amino acids in ATOM records.
                # Preserve coordinates and numbering; substitute MET for MSE.
                rid = int(s[22:26])
                if rid not in modifications:
                    modifications.append(rid)
                atom = ' SD ' if s[12:16].strip() == 'SE' else s[12:16]
                element = ' S' if s[12:16].strip() == 'SE' else s[76:78]
                s = 'ATOM  ' + s[6:12] + atom + s[16:17] + 'MET' + s[20:76] + element + s[78:]
            elif s.startswith('HETATM'):
                raise ValueError(f'Unsupported target HETATM: {pid} {s[17:20]}')
            atoms.append(s)
        coords, residues = {}, {}
        for s in atoms:
            xyz = tuple(float(s[i:i+8]) for i in (30, 38, 46))
            if not all(map(math.isfinite, xyz)) or s[26] != ' ':
                raise ValueError(f'Invalid target coordinates/numbering: {pid}')
            rid = int(s[22:26])
            residues.setdefault(rid, set()).add(s[12:16].strip())
            if s[12:16].strip() == 'CA':
                if rid in coords:
                    raise ValueError(f'Duplicate CA: {pid} {rid}')
                coords[rid] = xyz
        expected = {v for k, v in mapping.items() if k[0] == pid}
        if set(residues) != expected:
            raise ValueError(f'Target mapping mismatch: {pid}')
        candidates = [rid for _, rid in sorted(epitopes[pid])]
        if not set(candidates).issubset(coords):
            raise ValueError(f'Epitope lacks CA coordinates: {pid}')
        chosen = []
        for rid in candidates:
            if all(math.dist(coords[rid], coords[c]) >= 8 for c in chosen):
                chosen.append(rid)
            if len(chosen) == 4:
                break
        relaxed = len(chosen) < 4
        for rid in candidates:
            if len(chosen) >= 4:
                break
            if rid not in chosen:
                chosen.append(rid)
        missing = {str(k): sorted({'N', 'CA', 'C', 'O'} - v)
                   for k, v in residues.items() if not {'N', 'CA', 'C', 'O'} <= v}
        record = dict(pdb_id=pid, antigen_residues=len(coords),
                      hotspots=[f'T{x}' for x in sorted(chosen)],
                      epitope_T=candidates, spacing_relaxed=relaxed,
                      missing_backbone_atoms=missing, mse_to_met_residues=modifications, qc=qc[pid],
                      source_sha256=digest(source))
        prepared.append((record, '\n'.join(atoms) + '\nTER\nEND\n'))
    dest = PROC / 'rfab_inputs'
    (dest / 'targets').mkdir(parents=True, exist_ok=True)
    rows = []
    for record, content in sorted(prepared, key=lambda x: x[0]['antigen_residues']):
        target = dest / 'targets' / (record['pdb_id'] + '_T.pdb')
        target.write_text(content)
        record.update(target=str(target), target_sha256=digest(target))
        rows.append(record)
    write_json(dest / 'targets.json', rows)
    print(f'Prepared {len(rows)} targets; unmapped epitope residues: 0', flush=True)
    return rows


def validate(folder, stage, expected=2):
    files = sorted(folder.glob('*best.pdb' if stage == 'rf2' else '*.pdb'))
    if len(files) != expected:
        raise ValueError(f'{stage}: expected {expected} structures, found {len(files)}')
    info = [inspect_pdb(p) for p in files]
    if stage == 'rf2':
        for item in info:
            for name in ('pred_lddt', 'interaction_pae'):
                if item['scores'].get(name) is None:
                    raise ValueError(f'Missing confidence metric: {name}')
    return info


def report(run, states):
    rows = []
    for s in states:
        stages = s.get('stages', {})
        def count(stage):
            if stage in stages:
                return len(stages[stage]['outputs'])
            attempts = [a for a in s.get('attempts', []) if a['stage'] == stage]
            return len(attempts[-1].get('partial_outputs', [])) if attempts else 0
        rows.append(dict(pdb_id=s['pdb_id'], arm=s['arm'], status=s['status'],
            backbones=count('rfdiffusion'), sequences=count('proteinmpnn'), predictions=count('rf2'),
            wall_seconds=round(sum(a.get('seconds', 0) for a in s.get('attempts', [])), 2),
            error=s.get('error', '')))
    with (run / 'status.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    write_json(run / 'summary.json', dict(updated_utc=stamp(), runs=rows,
        wall_seconds_all_attempts=sum(r['wall_seconds'] for r in rows),
        complete=sum(r['status'] == 'complete' for r in rows), total=len(rows)))
    shutil.copyfile(run / 'status.csv', ROOT / 'results/generation_status.csv')
    candidates = []
    sequences = []
    aa = dict(pair.split(':') for pair in
              'ALA:A ARG:R ASN:N ASP:D CYS:C GLN:Q GLU:E GLY:G HIS:H ILE:I LEU:L LYS:K MET:M PHE:F PRO:P SER:S THR:T TRP:W TYR:Y VAL:V'.split())
    for s in states:
        for item in s.get('stages', {}).get('rf2', {}).get('outputs', []):
            candidates.append(dict(pdb_id=s['pdb_id'], arm=s['arm'], file=item['file'], **item['scores']))
        for item in s.get('stages', {}).get('proteinmpnn', {}).get('outputs', []):
            path = Path(item['file'])
            seq = ''.join(aa.get(line[17:20], 'X') for line in path.read_text().splitlines()
                          if line.startswith('ATOM') and line[21] == 'H' and line[12:16].strip() == 'CA')
            sequences.append(f'>{s["pdb_id"]}_{s["arm"]}_{path.stem}\n{seq}\n')
    (run / 'designed_sequences.fasta').write_text(''.join(sequences))
    write_json(run / 'candidate_metrics.json', candidates)
    if candidates:
        keys = sorted({key for row in candidates for key in row})
        with (run / 'candidate_metrics.csv').open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(candidates)
        shutil.copyfile(run / 'candidate_metrics.csv', ROOT / 'results/rfab_candidate_metrics.csv')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prepare-only', action='store_true')
    p.add_argument('--pilot-only', action='store_true')
    p.add_argument('--retry-failed', action='store_true')
    p.add_argument('--rfantibody-root', type=Path, default=ROOT.parent / 'RFantibody')
    args = p.parse_args()
    targets = prepare()
    if args.prepare_only:
        return
    up = args.rfantibody_root.resolve()
    env = dict(os.environ, PATH=str(up / '.venv/bin') + ':' + os.environ['PATH'],
               RFANTIBODY_ROOT=str(up), PYTHONUNBUFFERED='1', HYDRA_FULL_ERROR='1')
    subprocess.run(['python', '-c', 'import torch,dgl; assert torch.cuda.is_available(); '
                    'print(torch.__version__,dgl.__version__,torch.cuda.get_device_name(0))'],
                   env=env, check=True)
    framework = up / 'scripts/examples/example_inputs/h-NbBCII10.pdb'
    cfg = dict(version=1, loops='H1:7,H2:6,H3:8', seeds=[0, 1], temperature=0.1,
               rf2_recycles=1, rf2_seed=37, framework_sha256=digest(framework),
               upstream_commit=subprocess.check_output(['git', '-C', str(up), 'rev-parse', 'HEAD'], text=True).strip(),
               targets=[(r['pdb_id'], r['target_sha256'], r['hotspots']) for r in targets],
               weights={name: digest(up / 'weights' / name) for name in
                        ('RFdiffusion_Ab.pt', 'ProteinMPNN_v48_noise_0.2.pt', 'RF2_ab.pt')})
    runid = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:12]
    run = ROOT / 'runs' / ('rfab_benchmark_' + runid)
    run.mkdir(parents=True, exist_ok=True)
    import fcntl
    lock = (run / 'runner.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    write_json(run / 'config.json', cfg)
    write_json(ROOT / 'configs/rfab_benchmark.json', cfg)
    shutil.copytree(PROC / 'rfab_inputs', run / 'inputs', dirs_exist_ok=True)
    for script in ('benchmark_rfantibody.py', 'pilot_rfantibody.py'):
        shutil.copyfile(Path(__file__).parent / script, run / script)
    shutil.copyfile(framework, run / 'framework.pdb')
    framework = run / 'framework.pdb'
    states = []
    for target in targets:
        for arm in ('hotspot', 'unguided'):
            folder = run / (target['pdb_id'] + '_' + arm)
            folder.mkdir(exist_ok=True)
            statefile = folder / 'state.json'
            s = json.loads(statefile.read_text()) if statefile.exists() else dict(
                pdb_id=target['pdb_id'], arm=arm, status='unattempted', stages={}, attempts=[])
            states.append(s)
    report(run, states)
    print(f'Run directory: {run}', flush=True)
    for index, s in enumerate(states):
        if args.pilot_only and index > 0:
            break
        target = next(t for t in targets if t['pdb_id'] == s['pdb_id'])
        folder = run / (s['pdb_id'] + '_' + s['arm'])
        statefile = folder / 'state.json'
        if s['status'] == 'failed' and not args.retry_failed:
            continue
        try:
            s['status'] = 'running'
            s.pop('error', None)
            for stage in ('rfdiffusion', 'proteinmpnn', 'rf2'):
                if stage in s['stages']:
                    validate(Path(s['stages'][stage]['directory']), stage)
                    continue
                attempt = folder / f'{stage}_attempt_{sum(a["stage"] == stage for a in s["attempts"]) + 1}'
                attempt.mkdir()
                output = attempt / 'outputs'
                output.mkdir()
                if stage == 'rfdiffusion':
                    cmd = ['rfdiffusion', '-t', target['target'], '-f', str(framework),
                           '-o', str(output / 'nb'), '-n', '2', '-l', cfg['loops'],
                           '--diffuser-t', '50', '--deterministic', '--no-trajectory']
                    if s['arm'] == 'hotspot':
                        cmd += ['-h', ','.join(target['hotspots'])]
                elif stage == 'proteinmpnn':
                    cmd = ['proteinmpnn', '-i', s['stages']['rfdiffusion']['directory'],
                           '-o', str(output), '-n', '1', '-l', 'H1,H2,H3', '-t', '0.1',
                           '--deterministic', '--debug']
                else:
                    cmd = ['rf2', '-i', s['stages']['proteinmpnn']['directory'], '-o', str(output),
                           '-r', '1', '-s', '37', '--hotspot-show-prop', '0.0']
                a = dict(stage=stage, started_utc=stamp(), command=cmd, log=str(attempt / 'output.log'))
                s['attempts'].append(a)
                write_json(statefile, s)
                report(run, states)
                print(f'{index+1}/42 {s["pdb_id"]} {s["arm"]}: {stage}', flush=True)
                tick = time.monotonic()
                try:
                    with (attempt / 'output.log').open('w') as log:
                        proc = subprocess.Popen(cmd, env=env, cwd=attempt, stdout=log,
                                                stderr=subprocess.STDOUT, start_new_session=True)
                        try:
                            code = proc.wait(timeout=3600)
                        except BaseException:
                            os.killpg(proc.pid, signal.SIGTERM)
                            try:
                                proc.wait(timeout=10)
                            except subprocess.TimeoutExpired:
                                os.killpg(proc.pid, signal.SIGKILL)
                                proc.wait()
                            raise
                    a['returncode'] = code
                    if code:
                        raise RuntimeError(f'{stage} exit {code}; see {a["log"]}')
                    info = validate(output, stage)
                    s['stages'][stage] = dict(directory=str(output), outputs=info)
                finally:
                    a['partial_outputs'] = []
                    for path in sorted(output.glob('*best.pdb' if stage == 'rf2' else '*.pdb')):
                        try:
                            a['partial_outputs'].append(inspect_pdb(path))
                        except ValueError:
                            pass
                    a['seconds'] = round(time.monotonic() - tick, 2)
                    a['finished_utc'] = stamp()
                    write_json(statefile, s)
            s['status'] = 'complete'
            print(f'COMPLETE: {s["pdb_id"]} {s["arm"]}', flush=True)
        except KeyboardInterrupt:
            s['status'] = 'interrupted'
            s['error'] = 'Stopped by operator; resume with the same command.'
            print(f'INTERRUPTED: {s["pdb_id"]} {s["arm"]}; completed stages retained.', flush=True)
            raise SystemExit(130)
        except Exception as exc:
            s['status'] = 'failed'
            s['error'] = str(exc)
            print(f'FAILED: {s["pdb_id"]} {s["arm"]}: {exc}', flush=True)
        finally:
            write_json(statefile, s)
            report(run, states)
        if index == 0 and s['status'] != 'complete':
            raise SystemExit('Pilot failed; sweep not started. Inspect state.json and logs.')
    print(f'Results: {run / "summary.json"}', flush=True)


if __name__ == '__main__':
    main()
