#!/usr/bin/env python3
"""Validate AVIDa source tables; exact edit-distance components; whole-VHH splits."""
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
from rapidfuzz import process
from rapidfuzz.distance import Levenshtein

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src/data'))
from download_avida import digest

AA = set('ACDEFGHIKLMNPQRSTVWY')
ASSAY = 'released_phage_display_biopanning_NGS_binary_labels'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def stable_id(text):
    return hashlib.sha256(text.encode()).hexdigest()[:24]


def validate_sources(pairs, antigens):
    require(set(['VHH_sequence', 'Ag_label', 'label', 'subject_species', 'subject_name', 'subject_sex']) <= set(pairs), 'Missing source columns')
    require(not pairs.isna().any().any() and not antigens.isna().any().any(), 'Missing source values')
    require(set(pairs.label.unique()) <= {0, 1}, 'Nonbinary labels')
    require(not antigens.Ag_label.duplicated().any(), 'Duplicate antigen key')
    for seq in set(pairs.VHH_sequence) | set(antigens.Ag_sequence):
        require(bool(seq) and set(seq) <= AA, 'Noncanonical or empty sequence')
    require(set(pairs.Ag_label) <= set(antigens.Ag_label), 'Unmatched antigen label')
    key = ['VHH_sequence', 'Ag_label']
    require(pairs.groupby(key).label.nunique().max() == 1, 'Conflicting experimental pair labels')
    # A repeated biological pair with differing metadata cannot be silently merged.
    duplicate_rows = int(pairs.duplicated().sum())
    pairs = pairs.drop_duplicates().copy()
    require(not pairs.duplicated(key).any(), 'Repeated pair with differing subject/assay metadata')
    pairs['source_row'] = pairs.index + 2
    return pairs.merge(antigens, on='Ag_label', how='left', validate='many_to_one'), duplicate_rows


def components(sequences, block_size=256, workers=4, progress=True):
    """Exhaustive pairs: union iff 10*Levenshtein <= max length (exact 90% edit similarity)."""
    n = len(sequences)
    parent = np.arange(n)
    lengths = np.array(list(map(len, sequences)))
    cutoff = int(lengths.max() // 10)
    require(cutoff < 254, 'Sequence too long for bounded uint8 distances')
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return int(i)
    edges, comparisons = 0, 0
    for start in range(0, n, block_size):
        stop = min(n, start + block_size)
        distances = process.cdist(sequences[start:stop], sequences[start:], scorer=Levenshtein.distance,
                                 score_cutoff=cutoff, dtype=np.uint8, workers=workers)
        for offset, i in enumerate(range(start, stop)):
            js = np.flatnonzero(distances[offset].astype(np.int16) * 10 <= np.maximum(lengths[i], lengths[start:])) + start
            js = js[js > i]
            edges += len(js)
            comparisons += n - i - 1
            root = find(i)
            for j in js:
                other = find(j)
                if other != root:
                    parent[other] = root
        if progress and (start // block_size) % 10 == 0:
            print(f'Cluster comparisons: {stop}/{n} sequences', flush=True)
    roots = [find(i) for i in range(n)]
    names = {}
    for i, root in enumerate(roots):
        names[root] = min(names.get(root, sequences[i]), sequences[i])
    return ['vhh_cluster_' + stable_id(names[r]) for r in roots], dict(similarity_edges=edges, pairwise_comparisons=comparisons)


def split_components(vhhs, row_counts, seed=37):
    """Greedy row-balanced component allocation, using no labels; largest groups first."""
    sizes = vhhs.assign(rows=vhhs.vhh_id.map(row_counts)).groupby('cluster_id').rows.sum().to_dict()
    names, fractions = ['train', 'development', 'test'], [0.7, 0.15, 0.15]
    used = dict.fromkeys(names, 0)
    assignment = {}
    ordered = sorted(sizes, key=lambda c: (-sizes[c], stable_id(f'{seed}:{c}')))
    for cluster in ordered:
        dest = min(range(3), key=lambda j: used[names[j]] / fractions[j])
        assignment[cluster] = names[dest]
        used[names[dest]] += sizes[cluster]
    return vhhs.cluster_id.map(assignment)


def main():
    tick = time.monotonic()
    cfg = json.loads((ROOT / 'configs/avida_protocol.json').read_text())
    require(cfg['clustering']['maximum_distance_fraction'] == .1, 'This implementation fixes edit distance at 10%')
    source = json.loads((ROOT / 'results/evaluation/avida_source_manifest.json').read_text())
    require(source['revision'] == cfg['revision'], 'Revision mismatch')
    files = {r['name']: ROOT / r['path'] for r in source['files']}
    for row in source['files']:
        require(digest(ROOT / row['path']) == row['sha256'], 'Raw source checksum mismatch')
    raw = pd.read_csv(files['AVIDa-hIL6.csv'])
    pairs, duplicates = validate_sources(raw, pd.read_csv(files['antigen_sequences.csv']))
    pairs['dataset'] = 'AVIDa-hIL6'
    pairs['assay_group'] = ASSAY
    sequences = sorted(pairs.VHH_sequence.unique())
    vhhs = pd.DataFrame({'VHH_sequence': sequences})
    vhhs['vhh_id'] = ['vhh_' + stable_id(s) for s in sequences]
    require(vhhs.vhh_id.is_unique, 'VHH hash collision')
    vhhs['cluster_id'], cluster_stats = components(sequences, cfg['clustering']['block_size'], cfg['clustering']['workers'])
    mapping = dict(zip(vhhs.VHH_sequence, vhhs.vhh_id))
    pairs['vhh_id'] = pairs.VHH_sequence.map(mapping)
    vhhs['split'] = split_components(vhhs, pairs.vhh_id.value_counts(), cfg['split']['seed'])
    pairs = pairs.merge(vhhs[['vhh_id', 'cluster_id', 'split']], on='vhh_id', validate='many_to_one')
    pairs['experimental_pair_id'] = ['avida_' + stable_id(f'AVIDa-hIL6|{s}|{a}|{t}|{ASSAY}')
                                      for s, a, t in zip(pairs.VHH_sequence, pairs.Ag_label, pairs.Ag_sequence)]
    require(pairs.experimental_pair_id.is_unique, 'Pair identity collision')
    require(pairs.groupby('vhh_id').split.nunique().max() == 1, 'VHH split leakage')
    require(pairs.groupby('cluster_id').split.nunique().max() == 1, 'Cluster split leakage')
    require(len(pairs) == 573891 and pairs.label.sum() == 20980 and len(vhhs) == 38599 and pairs.Ag_label.nunique() == 31,
            'Unexpected release counts; investigate source rather than forcing counts')
    out = ROOT / 'data/processed'
    out.mkdir(exist_ok=True, parents=True)
    vhhs.to_csv(out / 'avida_vhh_splits.csv.gz', index=False, compression={'method':'gzip', 'mtime':0})
    pairs.to_csv(out / 'avida_pairs.csv.gz', index=False, compression={'method':'gzip', 'mtime':0})
    summary = []
    for split in ['all', 'train', 'development', 'test']:
        df = pairs if split == 'all' else pairs[pairs.split == split]
        for antigen in ['all'] + sorted(df.Ag_label.unique()):
            part = df if antigen == 'all' else df[df.Ag_label == antigen]
            summary.append(dict(split=split, antigen=antigen, pairs=len(part), positives=int(part.label.sum()),
                negatives=int((part.label == 0).sum()), prevalence=float(part.label.mean()),
                unique_vhhs=part.vhh_id.nunique(), clusters=part.cluster_id.nunique()))
    results = ROOT / 'results/evaluation'
    pd.DataFrame(summary).to_csv(results / 'avida_class_balance.csv', index=False)
    pilot_rows, used = [], set()
    pool = pairs[(pairs.split == cfg['pilot']['split']) & (pairs.Ag_label == cfg['pilot']['antigen'])].sort_values('experimental_pair_id')
    for label, number in [(1, cfg['pilot']['positives']), (0, cfg['pilot']['negatives'])]:
        selected = 0
        for _, row in pool[pool.label == label].iterrows():
            if row.cluster_id not in used:
                pilot_rows.append(row.to_dict())
                used.add(row.cluster_id)
                selected += 1
                if selected == number:
                    break
        require(selected == number, 'Insufficient distinct development clusters for balanced pilot')
    pd.DataFrame(pilot_rows).to_csv(results / 'avida_pilot_pairs.csv', index=False)
    audit = dict(source_rows=len(raw), duplicate_rows_removed=duplicates, curated_pairs=len(pairs),
        positives=int(pairs.label.sum()), negatives=int((pairs.label == 0).sum()),
        unique_vhhs=len(vhhs), clusters=vhhs.cluster_id.nunique(), antigens=pairs.Ag_label.nunique(),
        largest_cluster_vhhs=int(vhhs.cluster_id.value_counts().max()),
        largest_cluster_pairs=int(pairs.cluster_id.value_counts().max()),
        cross_split_vhh_overlap=0, cross_split_component_overlap=0,
        close_pair_cross_split_edges=0, **cluster_stats,
        seconds=round(time.monotonic()-tick, 3), protocol_sha256=digest(ROOT / 'configs/avida_protocol.json'),
        outputs={p.name:digest(p) for p in [out/'avida_vhh_splits.csv.gz', out/'avida_pairs.csv.gz']},
        metadata_limit='Released subject metadata retained; per-pair library, assay replicate, concentration and KD are unavailable',
        sequence_policy='Released full sequences retained verbatim, including leader-like prefixes and His6 tails; no undocumented trimming')
    (results / 'avida_curation_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps(audit, indent=2), flush=True)


if __name__ == '__main__':
    main()
