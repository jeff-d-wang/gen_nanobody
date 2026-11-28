#!/usr/bin/env python3
"""Download only ANDD v2 main XLSX and data dictionary; verify Zenodo checksums."""
from pathlib import Path
import hashlib, json, requests
from urllib.parse import urlparse

ROOT = Path(__file__).resolve().parent
RAW = ROOT / 'raw' / 'andd_v2'
RAW.mkdir(parents=True, exist_ok=True)
RECORD = 'https://zenodo.org/api/records/18151718'

def md5(path):
    h=hashlib.md5()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024), b''): h.update(chunk)
    return h.hexdigest()

def main():
    r=requests.get(RECORD,timeout=60);r.raise_for_status();record=r.json()
    (RAW/'zenodo_record.json').write_text(json.dumps(record,indent=2,ensure_ascii=False))
    files=record.get('files',[])
    selected=[]
    for f in files:
        key=f.get('key', '')
        if key.lower().endswith(('.xlsx','.csv')) and ('andd' in key.lower() or 'dictionary' in key.lower()):
            selected.append(f)
    if not selected or not any(x['key'].lower().endswith('.xlsx') for x in selected):
        raise RuntimeError('Expected ANDD XLSX not identified in Zenodo record; inspect raw/andd_v2/zenodo_record.json')
    manifest=[]
    for f in selected:
        key=f['key']; dest=RAW / Path(key).name
        url=(f.get('links') or {}).get('self') or f'https://zenodo.org/records/18151718/files/{key}?download=1'
        checksum=f.get('checksum','')
        if dest.exists() and (not checksum.startswith('md5:') or md5(dest)==checksum.split(':',1)[1]):
            print('CACHE',dest.name)
        else:
            print('DOWNLOAD',key)
            with requests.get(url,stream=True,timeout=(20,180)) as resp:
                resp.raise_for_status()
                with dest.with_suffix(dest.suffix+'.part').open('wb') as o:
                    for chunk in resp.iter_content(chunk_size=1024*1024):
                        if chunk:o.write(chunk)
            dest.with_suffix(dest.suffix+'.part').replace(dest)
        found=md5(dest)
        if checksum.startswith('md5:') and found!=checksum.split(':',1)[1]:
            raise RuntimeError(f'Checksum mismatch for {key}: {found} != {checksum}')
        manifest.append(dict(file=dest.name,zenodo_key=key,url=url,bytes=dest.stat().st_size,md5=found,expected_checksum=checksum))
    (RAW/'download_manifest.json').write_text(json.dumps(manifest,indent=2))
    print('Downloaded and checksum-validated',len(manifest),'files')

if __name__=='__main__':main()
