#!/usr/bin/env python3
"""Four exact released AVIDa pairs: RF2 throughput pilot, not an RQ1 comparison.

Run with the existing RFantibody Python. No upstream files are changed.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import urllib.request

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src/data'))
from download_avida import digest

AA = 'ARNDCQEGHILKMFPSTWYV'
AA3 = 'ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split()
TO1 = dict(zip(AA3, AA))
TEMPLATE_HASH = '37a666e3b62dca0dd61a029c430f936800efcb872f734655429816bfb9f845db'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def target_template(path, sequence):
    from rfantibody.rf2.network.chemical import aa2long
    lines = path.read_text().splitlines()
    dbrefs = [line for line in lines if line.startswith('DBREF ') and line[12] == 'A']
    require(len(dbrefs) == 1 and 'P05231' in dbrefs[0], 'Unexpected IL-6 DBREF')
    ref = dbrefs[0]
    pdb_start, pdb_end = int(ref[14:18]), int(ref[20:24])
    seq_start, seq_end = int(ref[55:60]), int(ref[62:67])
    require(seq_end - seq_start == pdb_end - pdb_start, 'Nonlinear DBREF mapping')
    offset = seq_start - pdb_start - 1
    xyz = np.full((len(sequence), 27, 3), np.nan, dtype=np.float32)
    mask = np.zeros((len(sequence), 27), dtype=bool)
    choices = {}
    for line in lines:
        if not line.startswith('ATOM') or line[21] != 'A':
            continue
        require(line[26] == ' ', 'Insertion code requires explicit target mapping')
        original = int(line[22:26])
        index = original + offset
        require(0 <= index < len(sequence), 'Template residue out of range')
        require(TO1[line[17:20]] == sequence[index], f'Target sequence mismatch at PDB residue {original}')
        atom = line[12:16].strip()
        names = [a.strip() if a else '' for a in aa2long[AA.index(sequence[index])]]
        if atom not in names:
            continue
        j = names.index(atom)
        occupancy = float(line[54:60])
        key = (index, j)
        if key not in choices or occupancy > choices[key]:
            choices[key] = occupancy
            xyz[index, j] = [float(line[k:k+8]) for k in (30, 38, 46)]
            mask[index, j] = True
    require(np.isfinite(xyz[mask]).all(), 'Nonfinite observed target coordinates')
    observed = np.flatnonzero(mask[:, :3].all(axis=1))
    require(len(observed) > 100, 'Insufficient target backbone coverage')
    return dict(seq=[AA.index(a) for a in sequence], xyz=xyz.tolist(), mask=mask.tolist(),
                idx=list(range(1, len(sequence)+1))), observed.tolist(), offset


def inspect_prediction(path, expected_h, expected_t, supported):
    atoms, sequences, scores = {'H': [], 'T': []}, {'H': [], 'T': []}, {}
    for line in path.read_text().splitlines():
        if line.startswith('ATOM'):
            chain = line[21]
            require(chain in atoms, f'Unexpected chain {chain}')
            xyz = [float(line[i:i+8]) for i in (30, 38, 46)]
            require(all(math.isfinite(v) for v in xyz), 'Nonfinite prediction coordinates')
            if line[12:16].strip() == 'CA':
                sequences[chain].append((line[22:27], TO1[line[17:20]]))
            if line[76:78].strip() != 'H' and not line[12:16].strip().startswith('H'):
                atoms[chain].append((line[22:27], xyz))
        elif line.startswith('SCORE '):
            key, value = line[6:].split(':', 1)
            scores[key] = float(value)
    for chain, expected in [('H', expected_h), ('T', expected_t)]:
        require(len(set(k for k, _ in sequences[chain])) == len(expected), 'Duplicate/missing residues')
        require(''.join(a for _, a in sequences[chain]) == expected, f'Prediction changed chain {chain}')
    require(all(math.isfinite(v) for v in scores.values()), 'Nonfinite confidence')
    require({'pred_lddt', 'interaction_pae'} <= set(scores), 'Missing confidence')
    require(0 <= scores['pred_lddt'] <= 1 and scores['interaction_pae'] >= 0, 'Invalid confidence range')
    t_indices = {key: i for i, (key, _) in enumerate(sequences['T'])}
    h_xyz = np.array([v for _, v in atoms['H']])
    contacts, severe = set(), 0
    t_xyz = np.array([v for _, v in atoms['T']])
    for start in range(0, len(h_xyz), 128):
        d2 = ((h_xyz[start:start+128, None, :] - t_xyz[None, :, :]) ** 2).sum(axis=2)
        contacts.update(t_indices[atoms['T'][j][0]] for j in np.flatnonzero((d2 <= 25).any(axis=0)))
        severe += int((d2 < 4).sum())
    return dict(**scores, target_contact_residues_5A=len(contacts),
                template_supported_target_contacts_5A=len(contacts & set(supported)),
                interchain_atom_pairs_below_2A=severe)


def write_rows(path, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rfantibody-root', type=Path, default=ROOT.parent / 'RFantibody')
    parser.add_argument('--out', type=Path, default=ROOT / 'runs/avida_pilot')
    args = parser.parse_args()
    up, out = args.rfantibody_root.resolve(), args.out.resolve()
    require(not out.exists() or not any(out.iterdir()), 'Use a new/empty pilot output directory')
    out.mkdir(exist_ok=True, parents=True)
    results = ROOT / 'results/evaluation'
    with (results / 'avida_pilot_pairs.csv').open() as stream:
        pairs = list(csv.DictReader(stream))
    require(len(pairs) == 4 and sorted(int(r['label']) for r in pairs) == [0, 0, 1, 1], 'Expected four balanced pairs')
    require(all(r['split'] == 'development' and r['Ag_label'] == 'IL-6_WTs' for r in pairs), 'Pilot must use development WT only')
    require(len({r['cluster_id'] for r in pairs}) == 4, 'Pilot clusters must differ')
    source = ROOT / 'data/raw/avida_pilot/1ALU.pdb'
    source.parent.mkdir(parents=True, exist_ok=True)
    if not source.exists():
        source.write_bytes(urllib.request.urlopen('https://files.rcsb.org/download/1ALU.pdb', timeout=30).read())
    require(digest(source) == TEMPLATE_HASH, 'Public target structure hash changed; review before use')
    target, supported, offset = target_template(source, pairs[0]['Ag_sequence'])
    import torch
    preflight = dict(cuda_available=torch.cuda.is_available(), torch_version=torch.__version__,
        gpu=torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
        gpu_memory_free_total_bytes=list(torch.cuda.mem_get_info()) if torch.cuda.is_available() else None,
        boltz_installed=importlib.util.find_spec('boltz') is not None,
        pyrosetta_installed=importlib.util.find_spec('pyrosetta') is not None,
        upstream_commit=subprocess.check_output(['git', '-C', str(up), 'rev-parse', 'HEAD'], text=True).strip(),
        upstream_dirty=bool(subprocess.check_output(['git', '-C', str(up), 'status', '--porcelain'], text=True)),
        rf2_weights_sha256=digest(up / 'weights/RF2_ab.pt'),
        pilot_pairs_sha256=digest(results / 'avida_pilot_pairs.csv'),
        template_url='https://files.rcsb.org/download/1ALU.pdb', template_sha256=TEMPLATE_HASH,
        template_backbone_residues=len(supported), full_released_target_residues=len(pairs[0]['Ag_sequence']),
        pdb_to_release_zero_based_offset=offset, missing_target_residues=len(pairs[0]['Ag_sequence'])-len(supported),
        template_policy='Observed chain A residues must match exact released WT sequence via DBREF; all missing atoms masked; full released sequences retained',
        parameters=dict(recycles=1, seed=37, hotspot_show_proportion=0., timeout_seconds_per_pair=600),
        limitations='Pilot only; RF2 not independent final predictor; signal/His-tag construct policy unvalidated; no Rosetta score; no balanced-pilot AUPRC')
    (out / 'preflight.json').write_text(json.dumps(preflight, indent=2) + '\n')
    (results / 'avida_pilot_preflight.json').write_text(json.dumps(preflight, indent=2) + '\n')
    env = dict(os.environ, PATH=str(up / '.venv/bin') + os.pathsep + os.environ['PATH'],
               RFANTIBODY_ROOT=str(up), PYTHONUNBUFFERED='1')
    rows = []
    for pair in pairs:
        identity = pair['experimental_pair_id']
        folder = out / identity
        folder.mkdir()
        inputs = folder / 'input.json'
        # Python/upstream RF2 JSON supports NaN placeholders for explicitly masked atoms.
        inputs.write_text(json.dumps({identity:dict(Hseq=pair['VHH_sequence'], Lseq='',
                              hotspots=[False]*len(pair['Ag_sequence']), T=target)}))
        output = folder / 'outputs'
        command = [str(up / '.venv/bin/rf2'), '-j', str(inputs), '-o', str(output),
                   '-r', '1', '-s', '37', '--hotspot-show-prop', '0.0']
        row = dict(experimental_pair_id=identity, cluster_id=pair['cluster_id'], split=pair['split'],
            antigen=pair['Ag_label'], label=int(pair['label']), status='not_started', model_process_wall_seconds='',
            structure='', error='', confidence_available=False, geometry_available=False,
            energy_available=False, independent_final_prediction=False)
        tick = time.monotonic()
        try:
            require(preflight['cuda_available'], 'CUDA unavailable')
            with (folder / 'model.log').open('w') as log:
                process = subprocess.Popen(command, env=env, cwd=folder, stdout=log, stderr=subprocess.STDOUT,
                                           start_new_session=True)
                try:
                    code = process.wait(timeout=600)
                except BaseException:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                    raise
            row['model_process_wall_seconds'] = round(time.monotonic()-tick, 3)
            require(code == 0, f'RF2 exit {code}; see model.log')
            paths = list(output.glob('*best.pdb'))
            require(len(paths) == 1, 'Expected one best prediction')
            metrics = inspect_prediction(paths[0], pair['VHH_sequence'], pair['Ag_sequence'], supported)
            row.update(status='passed', structure=paths[0].relative_to(ROOT).as_posix(),
                       confidence_available=True, geometry_available=True, **metrics)
        except Exception as exc:
            row.update(status='failed', error=f'{type(exc).__name__}: {exc}')
            if row['model_process_wall_seconds'] == '':
                row['model_process_wall_seconds'] = round(time.monotonic()-tick, 3)
        rows.append(row)
        write_rows(results / 'avida_feasibility.csv', rows)
        (folder / 'attempt.json').write_text(json.dumps(dict(result=row, command=command,
            input_sha256=digest(inputs)), indent=2) + '\n')
        print(json.dumps(row), flush=True)
    summary = dict(pilot_pairs=len(rows), passed=sum(r['status']=='passed' for r in rows),
        total_model_process_wall_seconds=round(sum(r['model_process_wall_seconds'] for r in rows), 3),
        full_structure_rq1_gate='no_go',
        reason='Four balanced development pairs are only a throughput pilot. Energy and an independent final predictor are unavailable; released construct handling needs validation before a structural benchmark.',
        no_large_modeling_campaign_started=True)
    (results / 'avida_feasibility_decision.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
