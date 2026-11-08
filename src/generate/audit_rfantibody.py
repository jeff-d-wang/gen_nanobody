#!/usr/bin/env python3
"""Audit a completed benchmark and write a portable handoff summary.

Usage: python3 src/generate/audit_rfantibody.py runs/rfab_benchmark_<hash>
Does not run models or change their outputs. Exits nonzero if anything is
unattempted, incomplete, failed, or inconsistent with the recorded inputs.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

from benchmark_rfantibody import ROOT, validate, write_json


def residues(path, chain):
    return [s[17:20] for s in path.read_text().splitlines()
            if s.startswith('ATOM') and s[21] == chain and s[12:16].strip() == 'CA']


def check_lineage(backbone, designed, predicted, target):
    """Check antigen identity, fixed framework, and sequence preservation."""
    target_seq = residues(target, 'T')
    for path in (backbone, designed, predicted):
        if residues(path, 'T') != target_seq:
            raise ValueError(f'Target sequence changed: {path}')
    before, after = residues(backbone, 'H'), residues(designed, 'H')
    if len(before) != len(after):
        raise ValueError(f'Heavy-chain length changed: {designed}')
    cdr = set()
    for line in backbone.read_text().splitlines():
        if line.startswith('REMARK PDBinfo-LABEL:'):
            fields = line.split()
            if fields[3] in ('H1', 'H2', 'H3'):
                cdr.add(int(fields[2]) - 1)
    if not cdr or not cdr <= set(range(len(before))):
        raise ValueError(f'Invalid CDR indices: {backbone}')
    if any(a != b for i, (a, b) in enumerate(zip(before, after)) if i not in cdr):
        raise ValueError(f'ProteinMPNN changed the framework: {designed}')
    if residues(predicted, 'H') != after:
        raise ValueError(f'RF2 changed the designed sequence: {predicted}')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('run', type=Path)
    args = p.parse_args()
    run = args.run.resolve()
    config = json.loads((run / 'config.json').read_text())
    rows, errors = [], []
    for pid, target_hash, _ in config['targets']:
        target = run / 'inputs/targets' / f'{pid}_T.pdb'
        if hashlib.sha256(target.read_bytes()).hexdigest() != target_hash:
            errors.append(f'Input hash mismatch: {pid}')
        for arm in ('hotspot', 'unguided'):
            folder = run / f'{pid}_{arm}'
            statefile = folder / 'state.json'
            state = json.loads(statefile.read_text()) if statefile.exists() else {}
            row = dict(pdb_id=pid, arm=arm, status=state.get('status', 'unattempted'),
                       validated_candidates=0, error=state.get('error', ''))
            if row['status'] == 'complete':
                try:
                    outputs = {}
                    for stage, info in state['stages'].items():
                        # Resolve by directory name so a copied archive remains auditable.
                        output = folder / Path(info['directory']).parent.name / 'outputs'
                        outputs[stage] = validate(output, stage)
                    for item in outputs['rfdiffusion']:
                        bb = Path(item['file'])
                        seqdir = Path(outputs['proteinmpnn'][0]['file']).parent
                        rfdir = Path(outputs['rf2'][0]['file']).parent
                        seq = seqdir / f'{bb.stem}_dldesign_0.pdb'
                        rf = rfdir / f'{bb.stem}_dldesign_0_best.pdb'
                        check_lineage(bb, seq, rf, target)
                        row['validated_candidates'] += 1
                except Exception as exc:
                    row.update(status='audit_failed', error=str(exc))
            if row['status'] != 'complete':
                errors.append(f'{pid} {arm}: {row["status"]} {row["error"]}')
            rows.append(row)
    audit = dict(passed=not errors, target_conditions=len(rows),
                 validated_candidates=sum(r['validated_candidates'] for r in rows),
                 errors=errors, rows=rows)
    write_json(run / 'audit.json', audit)
    with (run / 'audit.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = json.loads((run / 'summary.json').read_text())
    complete = sum(r['status'] == 'complete' for r in rows)
    text = (f'# RFantibody benchmark execution\n\n'
            f'- Run: `{run.name}`\n'
            f'- Complete conditions: {complete}/{len(rows)} (21 targets, two conditions).\n'
            f'- Candidates passing lineage checks: {audit["validated_candidates"]}/84.\n'
            f'- Recorded model-process time, including failed attempts: '
            f'{summary["wall_seconds_all_attempts"] / 3600:.2f} hours.\n\n'
            'Settings: two deterministic backbones per target/condition; CDR lengths '
            '7/6/8; 50 diffusion steps; one ProteinMPNN sequence per backbone at '
            'temperature 0.1; RF2 with one recycle and seed 37.\n\n'
            'Read `status.csv` for progress and failures, `candidate_metrics.csv` '
            'for RF2 scores, `designed_sequences.fasta` for sequences, and '
            '`audit.json` for validation details. Each target/condition folder '
            'contains per-stage attempt logs, structures, and `state.json`. '
            '`config.json`, `inputs/`, and `framework.pdb` identify the inputs.\n\n'
            'Both conditions use the same reference-derived antigen crop. Only '
            'the hotspot condition receives explicit hotspot labels. This is a '
            'workflow feasibility run; no binding claims or threshold calibration '
            'are supported. RF2 is a pipeline predictor, not independent final scoring.\n\n'
            'Raw confidence matrices are not exported by this runner; the saved '
            'RF2 PDB score lines and tables contain aggregate confidence metrics.\n')
    if errors:
        text += '\nOutstanding issues:\n\n' + '\n'.join('- ' + e for e in errors) + '\n'
    (run / 'HANDOFF.md').write_text(text)
    (ROOT / 'results/RFANTIBODY_RUN.md').write_text(text)
    print(json.dumps({k: v for k, v in audit.items() if k != 'rows'}, indent=2))
    return 0 if audit['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
