#!/usr/bin/env python3
"""Audit the frozen RFantibody archive and report recorded model-process costs.

Standard library only. No models are run and no archive files are modified.
"""
import argparse
import csv
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src/generate'))
from audit_rfantibody import check_lineage, residues
from benchmark_rfantibody import validate

STAGES = ('rfdiffusion', 'proteinmpnn', 'rf2')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def portable(path):
    return PurePosixPath(str(path).replace('\\', '/'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def attempt_rows(state, accepted):
    rows = []
    names = set()
    for attempt in state['attempts']:
        stage = attempt['stage']
        require(stage in STAGES, f'Unknown stage: {stage}')
        name = portable(attempt['log']).parent.name
        require(name not in names, f'Duplicate attempt: {name}')
        names.add(name)
        seconds = float(attempt['seconds'])
        require(math.isfinite(seconds) and seconds >= 0, f'Invalid cost: {name}')
        start = datetime.fromisoformat(attempt['started_utc'])
        finish = datetime.fromisoformat(attempt['finished_utc'])
        require(start.tzinfo is not None and finish.tzinfo is not None and finish >= start,
                f'Invalid timestamps: {name}')
        is_accepted = name == accepted[stage]
        require(not is_accepted or attempt.get('returncode') == 0,
                f'Accepted attempt did not exit successfully: {name}')
        rows.append(dict(pdb_id=state['pdb_id'], arm=state['arm'], stage=stage,
                         attempt=name, accepted=is_accepted,
                         returncode=attempt.get('returncode', ''),
                         started_utc=attempt['started_utc'], finished_utc=attempt['finished_utc'],
                         wall_seconds=seconds, audited_candidates_added=2 if is_accepted and stage == 'rf2' else 0))
    require(all(sum(r['accepted'] and r['stage'] == s for r in rows) == 1 for s in STAGES),
            'Each stage must have exactly one accepted attempt')
    return rows


def analyze(run):
    config = json.loads((run / 'config.json').read_text())
    status = read_csv(ROOT / 'results/generation_status.csv')
    status_map = {(r['pdb_id'], r['arm']): r for r in status}
    require(len(status_map) == len(status) == 42, 'Expected 42 unique status rows')
    conditions, attempts, candidates, inputs = [], [], [], []
    sequences = {}
    def record(path):
        inputs.append(dict(path=path.relative_to(ROOT).as_posix(), sha256=sha(path)))
    record(run / 'config.json')
    record(ROOT / 'results/generation_status.csv')
    record(run / 'summary.json')
    require(len(config['targets']) == 21 and len({t[0] for t in config['targets']}) == 21,
            'Expected 21 distinct targets')
    for pid, target_hash, _ in config['targets']:
        target = run / 'inputs/targets' / f'{pid}_T.pdb'
        require(sha(target) == target_hash, f'Target hash mismatch: {pid}')
        record(target)
        for arm in ('hotspot', 'unguided'):
            folder = run / f'{pid}_{arm}'
            statefile = folder / 'state.json'
            record(statefile)
            state = json.loads(statefile.read_text())
            require((state['pdb_id'], state['arm'], state['status']) == (pid, arm, 'complete'),
                    f'Incomplete or mismatched condition: {folder.name}')
            require(set(state['stages']) == set(STAGES), 'Unexpected stage set')
            accepted, outputs = {}, {}
            for stage in STAGES:
                info = state['stages'][stage]
                name = portable(info['directory']).parent.name
                require(name.startswith(stage + '_attempt_'), f'Invalid accepted directory: {name}')
                accepted[stage] = name
                output = folder / name / 'outputs'
                outputs[stage] = [Path(i['file']) for i in validate(output, stage)]
                require({p.name for p in outputs[stage]} ==
                        {portable(i['file']).name for i in info['outputs']}, 'State/output mismatch')
                for path in outputs[stage]:
                    record(path)
            local_attempts = attempt_rows(state, accepted)
            attempts.extend(local_attempts)
            local_sequences = set()
            for backbone in outputs['rfdiffusion']:
                index = int(backbone.stem.removeprefix('nb_'))
                designed = folder / accepted['proteinmpnn'] / 'outputs' / f'{backbone.stem}_dldesign_0.pdb'
                predicted = folder / accepted['rf2'] / 'outputs' / f'{backbone.stem}_dldesign_0_best.pdb'
                require(designed in outputs['proteinmpnn'] and predicted in outputs['rf2'], 'Broken lineage')
                check_lineage(backbone, designed, predicted, target)
                # Exact full heavy-chain sequence; triplet residue names avoid lossy X conversion.
                sequence = '-'.join(residues(designed, 'H'))
                local_sequences.add(sequence)
                key = (pid, arm, index)
                require(key not in sequences, 'Duplicate candidate identity')
                sequences[key] = sequence
                candidates.append(dict(pdb_id=pid, arm=arm, backbone_index=index,
                    heavy_sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                    backbone=backbone.relative_to(ROOT).as_posix(),
                    designed=designed.relative_to(ROOT).as_posix(),
                    prediction=predicted.relative_to(ROOT).as_posix(), lineage_audit_passed=True))
            cost = round(sum(r['wall_seconds'] for r in local_attempts), 2)
            old = status_map[(pid, arm)]
            require(old['status'] == 'complete' and all(int(old[k]) == 2 for k in ('backbones', 'sequences', 'predictions')),
                    f'Status count mismatch: {pid} {arm}')
            require(abs(cost - float(old['wall_seconds'])) < .005, 'Status timing mismatch')
            conditions.append(dict(pdb_id=pid, arm=arm, completed_conditions=1,
                generated_backbones=2, designed_sequence_records=2, unique_heavy_sequences=len(local_sequences),
                audited_rf2_predictions=2, model_process_wall_seconds=cost,
                unaccepted_attempt_wall_seconds=round(sum(r['wall_seconds'] for r in local_attempts if not r['accepted']), 2)))
    require(len(candidates) == 84, 'Expected 84 audited candidates')
    total = round(sum(r['wall_seconds'] for r in attempts), 2)
    require(abs(total - json.loads((run / 'summary.json').read_text())['wall_seconds_all_attempts']) < .005,
            'Archived summary timing mismatch')
    summaries = []
    for scope in ('all', 'hotspot', 'unguided'):
        selected = [r for r in conditions if scope == 'all' or r['arm'] == scope]
        costs = [r for r in attempts if scope == 'all' or r['arm'] == scope]
        seconds = round(sum(r['wall_seconds'] for r in costs), 2)
        row = dict(scope=scope, completed_conditions=len(selected),
            generated_backbones=sum(r['generated_backbones'] for r in selected),
            designed_sequence_records=sum(r['designed_sequence_records'] for r in selected),
            unique_heavy_sequences=len({s for (_, arm, _), s in sequences.items() if scope == 'all' or arm == scope}),
            audited_rf2_predictions=sum(r['audited_rf2_predictions'] for r in selected),
            computational_filter_pass_candidates='', computational_filter_status='not_defined',
            model_process_wall_seconds=seconds, model_process_wall_hours=round(seconds / 3600, 8),
            attempts=len(costs), unaccepted_attempts=sum(not r['accepted'] for r in costs),
            unaccepted_attempt_wall_seconds=round(sum(r['wall_seconds'] for r in costs if not r['accepted']), 2))
        for stage in STAGES:
            row[stage + '_wall_seconds'] = round(sum(r['wall_seconds'] for r in costs if r['stage'] == stage), 2)
        summaries.append(row)
    attempts.sort(key=lambda r: (datetime.fromisoformat(r['finished_utc']), r['pdb_id'], r['arm'], r['attempt']))
    cumulative_seconds, cumulative_candidates = 0., 0
    for row in attempts:
        cumulative_seconds += row['wall_seconds']
        cumulative_candidates += row['audited_candidates_added']
        row['cumulative_wall_seconds'] = round(cumulative_seconds, 2)
        row['cumulative_audited_candidates'] = cumulative_candidates
    return summaries, conditions, attempts, candidates, inputs


def plot(path, attempts):
    """Dependency-free standalone SVG; coordinates are derived from saved attempt rows."""
    hours = attempts[-1]['cumulative_wall_seconds'] / 3600
    count = attempts[-1]['cumulative_audited_candidates']
    points = ['80,360']
    previous = 0
    for row in attempts:
        x = 80 + 640 * row['cumulative_wall_seconds'] / 3600 / hours
        y = 360 - 280 * row['cumulative_audited_candidates'] / count
        points.extend([f'{x:.2f},{360 - 280 * previous / count:.2f}', f'{x:.2f},{y:.2f}'])
        previous = row['cumulative_audited_candidates']
    ticks = ''.join(f'<text x="{80+640*i/4}" y="382" text-anchor="middle">{hours*i/4:.2f}</text>' for i in range(5))
    ticks += ''.join(f'<text x="68" y="{365-280*i/4}" text-anchor="end">{int(count*i/4)}</text>' for i in range(5))
    path.write_text(f'''<svg xmlns="http://www.w3.org/2000/svg" width="800" height="460" viewBox="0 0 800 460">
<rect width="800" height="460" fill="white"/>
<g font-family="Arial, sans-serif" font-size="14" fill="#182b3a">
<text x="400" y="28" text-anchor="middle" font-size="21">RFantibody execution yield and recorded cost</text>
<text x="400" y="52" text-anchor="middle">All attempts included; candidates counted after accepted RF2 completion</text>
<path d="M80 80 V360 H720" fill="none" stroke="#64748b"/>{ticks}
<polyline points="{' '.join(points)}" fill="none" stroke="#147d92" stroke-width="2.5"/>
<text x="400" y="414" text-anchor="middle">Cumulative model-process wall time (hours)</text>
<text transform="translate(22 220) rotate(-90)" text-anchor="middle">Audited RF2 candidates</text>
<text x="400" y="444" text-anchor="middle">Not GPU-hours or binding yield; ordered by recorded attempt completion</text>
</g></svg>''')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, default=ROOT / 'runs/rfab_benchmark_1ce3344e4822')
    parser.add_argument('--out', type=Path, default=ROOT / 'results/evaluation')
    args = parser.parse_args()
    summaries, conditions, attempts, candidates, inputs = analyze(args.run.resolve())
    args.out.mkdir(parents=True, exist_ok=True)
    for name, rows in [('summary', summaries), ('conditions', conditions), ('attempts', attempts), ('candidates', candidates), ('inputs', inputs)]:
        write_csv(args.out / f'rfab_compute_{name}.csv', rows)
    plot(args.out / 'rfab_compute_yield.svg', attempts)
    provenance = dict(source_commit=subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip(),
        script_sha256=sha(Path(__file__)), archive=args.run.resolve().relative_to(ROOT).as_posix(),
        input_manifest_sha256=sha(args.out / 'rfab_compute_inputs.csv'),
        audit='Revalidated accepted structures and sequence lineage without modifying archive',
        sequence_deduplication='Exact full H-chain residue sequence; all-scope deduplication is global',
        timing='Recorded seconds summed across all stage attempts, including unaccepted attempts',
        scientific_filter='Not defined; pass count missing, not zero')
    (args.out / 'rfab_compute_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps(summaries, indent=2))


if __name__ == '__main__':
    main()
