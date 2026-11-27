#!/usr/bin/env python3
"""Download the official AVIDa release at the pinned revision; verify on reuse."""
import hashlib
import json
from pathlib import Path
import shutil
import urllib.request

ROOT = Path(__file__).resolve().parents[2]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    protocol = json.loads((ROOT / 'configs/avida_protocol.json').read_text())
    out = ROOT / 'data/raw/avida_hil6' / protocol['revision']
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = ROOT / 'results/evaluation/avida_source_manifest.json'
    previous = json.loads(manifest_path.read_text()) if manifest_path.exists() else None
    if previous and (previous['dataset'], previous['revision']) != (protocol['dataset'], protocol['revision']):
        raise ValueError('Existing source manifest is for another release; review version changes explicitly')
    records = []
    for name in ('AVIDa-hIL6.csv', 'antigen_sequences.csv', 'README.md'):
        url = f'https://huggingface.co/datasets/{protocol["dataset"]}/resolve/{protocol["revision"]}/{name}'
        path = out / name
        expected = next((r['sha256'] for r in previous['files'] if r['name'] == name), None) if previous else None
        if not path.exists():
            print(f'Downloading {name}', flush=True)
            temp = path.with_suffix(path.suffix + '.part')
            with urllib.request.urlopen(url, timeout=120) as response, temp.open('wb') as stream:
                shutil.copyfileobj(response, stream)
            if expected and digest(temp) != expected:
                raise ValueError(f'Download hash mismatch: {name}')
            temp.replace(path)
        checksum = digest(path)
        if expected and checksum != expected:
            raise ValueError(f'Cached source hash mismatch: {name}')
        records.append(dict(name=name, url=url, path=path.relative_to(ROOT).as_posix(),
                            bytes=path.stat().st_size, sha256=checksum))
        print(f'{name}: {path.stat().st_size} bytes; SHA256 {checksum}', flush=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(dict(dataset=protocol['dataset'], revision=protocol['revision'],
        license=protocol['license'], paper='https://doi.org/10.48550/arXiv.2306.03329', files=records), indent=2) + '\n')


if __name__ == '__main__':
    main()
