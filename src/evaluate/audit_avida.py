#!/usr/bin/env python3
"""Create or verify the A2 artifact manifest; does not rerun models or fit classifiers."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src/data'))
from download_avida import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-manifest', action='store_true', help='Record current reviewed artifacts instead of verifying the saved manifest')
    args = parser.parse_args()
    destination = ROOT / 'results/evaluation/avida_artifact_manifest.json'
    if args.write_manifest:
        paths = [ROOT / 'A2_AVIDA_IMPLEMENTATION_REPORT.md', ROOT / 'configs/avida_protocol.json',
                 ROOT / 'src/data/download_avida.py', ROOT / 'src/data/prepare_avida.py', Path(__file__)]
        for pattern in ['avida_*.py', 'test_avida*.py', 'requirements-avida*.txt']:
            paths.extend((ROOT / 'src/evaluate').glob(pattern))
        paths.extend(p for p in destination.parent.glob('avida_*') if p != destination and p.is_file())
        paths.extend((ROOT / 'data/processed').glob('avida_*'))
        source = json.loads((destination.parent / 'avida_source_manifest.json').read_text())
        paths.extend(ROOT / row['path'] for row in source['files'])
        paths.append(ROOT / 'data/raw/avida_pilot/1ALU.pdb')
        for directory in ['runs/avida_baseline', 'runs/avida_pilot']:
            paths.extend(p for p in (ROOT / directory).rglob('*') if p.is_file())
        rows = [dict(path=p.relative_to(ROOT).as_posix(), bytes=p.stat().st_size, sha256=digest(p))
                for p in sorted(set(p.resolve() for p in paths))]
        content = dict(source_commit=subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip(),
                       report_date='2026-11-26', files=rows,
                       note='Raw data, curated full tables, models and pilot structures are ignored by Git; regenerate or transfer separately.')
        destination.write_text(json.dumps(content, indent=2) + '\n')
    data = json.loads(destination.read_text())
    errors = []
    for row in data['files']:
        path = ROOT / row['path']
        if not path.is_file() or path.stat().st_size != row['bytes'] or digest(path) != row['sha256']:
            errors.append(row['path'])
    print(json.dumps(dict(passed=not errors, files_checked=len(data['files']), errors=errors), indent=2))
    return bool(errors)


if __name__ == '__main__':
    raise SystemExit(main())
