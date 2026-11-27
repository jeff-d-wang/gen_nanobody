#!/usr/bin/env python3
"""Prespecified sequence-only AVIDa baseline, natural prevalence, cluster bootstrap."""
import gc
import json
from pathlib import Path
import sys
import time
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
import sklearn
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import average_precision_score, roc_auc_score, precision_recall_curve
from sklearn.preprocessing import OneHotEncoder

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src/data'))
from download_avida import digest
from prepare_avida import require


def metric_values(y, scores, k, weights=None):
    w = np.ones(len(y), dtype=np.int64) if weights is None else np.asarray(weights)
    active = w > 0
    require(len(np.unique(y[active])) == 2, 'Both classes required for discrimination metrics')
    # Input rows have pair-ID order; stable sorting gives label-independent tie resolution.
    order = np.argsort(-scores, kind='stable')
    before = np.cumsum(w[order]) - w[order]
    used = np.minimum(w[order], np.maximum(k - before, 0))
    denom = min(k, w.sum())
    return dict(average_precision=float(average_precision_score(y, scores, sample_weight=w)),
                auroc=float(roc_auc_score(y, scores, sample_weight=w)),
                precision_at_k=float(np.dot(used, y[order]) / denom),
                prevalence=float(np.dot(w, y) / w.sum()))


def bootstrap(y, scores, clusters, k, repetitions, seed):
    unique, inverse = np.unique(clusters, return_inverse=True)
    rng = np.random.default_rng(seed)
    samples = []
    for _ in range(repetitions):
        counts = np.bincount(rng.integers(0, len(unique), size=len(unique)), minlength=len(unique))
        weights = counts[inverse]
        if len(np.unique(y[weights > 0])) < 2:
            continue
        samples.append(metric_values(y, scores, k, weights))
    require(len(samples) >= repetitions * .95, 'Too many one-class bootstrap draws')
    return {key: np.quantile([s[key] for s in samples], [.025, .975]).tolist() for key in samples[0]}, len(samples)


def main():
    started = time.monotonic()
    cfg_path = ROOT / 'configs/avida_protocol.json'
    cfg = json.loads(cfg_path.read_text())
    out = ROOT / 'results/evaluation'
    audit = json.loads((out / 'avida_curation_audit.json').read_text())
    require(digest(cfg_path) == audit['protocol_sha256'], 'Protocol changed after curation; recurate explicitly')
    for filename, checksum in audit['outputs'].items():
        require(digest(ROOT / 'data/processed' / filename) == checksum, 'Curated input hash mismatch')
    pairs = pd.read_csv(ROOT / 'data/processed/avida_pairs.csv.gz').sort_values('experimental_pair_id').reset_index(drop=True)
    vhhs = pd.read_csv(ROOT / 'data/processed/avida_vhh_splits.csv.gz')
    train = pairs[pairs.split == 'train']
    require(set(train.label) == {0, 1}, 'Train requires both classes')
    vectorizer = TfidfVectorizer(analyzer='char', ngram_range=(3, 3), lowercase=False, dtype=np.float32)
    vectorizer.fit(vhhs.loc[vhhs.split == 'train', 'VHH_sequence'])
    vhh_features = vectorizer.transform(vhhs.VHH_sequence)
    index = dict(zip(vhhs.vhh_id, range(len(vhhs))))
    encoder = OneHotEncoder(handle_unknown='ignore', dtype=np.float32)
    encoder.fit(train[['Ag_label']])
    def features(frame):
        return sparse.hstack([vhh_features[frame.vhh_id.map(index).to_numpy()],
                              encoder.transform(frame[['Ag_label']])], format='csr', dtype=np.float32)
    b = cfg['baseline']
    model = SGDClassifier(loss=b['loss'], alpha=b['alpha'], max_iter=b['max_iter'], tol=b['tol'],
                          random_state=b['seed'], class_weight=b['class_weight'])
    tick = time.monotonic()
    x = features(train)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        model.fit(x, train.label.to_numpy())
    fit_seconds = time.monotonic() - tick
    del x
    gc.collect()
    artifacts = ROOT / 'runs/avida_baseline'
    artifacts.mkdir(exist_ok=True, parents=True)
    # Data-only model serialization: no pickle execution required to inspect coefficients.
    np.savez_compressed(artifacts / 'model.npz', coef=model.coef_, intercept=model.intercept_, idf=vectorizer.idf_)
    (artifacts / 'features.json').write_text(json.dumps(dict(vocabulary=vectorizer.vocabulary_, antigen_categories=encoder.categories_[0].tolist()), default=int))
    rows, predictions = [], []
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.3), sharex=True, sharey=True)
    e = cfg['evaluation']
    for ax, split in zip(axes, ['development', 'test']):
        frame = pairs[pairs.split == split].copy()
        tick = time.monotonic()
        x = features(frame)
        scores = model.predict_proba(x)[:, 1]
        score_seconds = time.monotonic() - tick
        del x
        y = frame.label.to_numpy()
        require(np.isfinite(scores).all() and ((scores >= 0) & (scores <= 1)).all(), 'Invalid scores')
        values = metric_values(y, scores, e['precision_k'])
        intervals, draws = bootstrap(y, scores, frame.cluster_id.to_numpy(), e['precision_k'],
                                    e['bootstrap_replicates'], e['bootstrap_seed'])
        for metric, value in values.items():
            rows.append(dict(model='VHH_trigram_plus_antigen_additive_logistic', split=split,
                metric=metric, value=value, ci95_low=intervals[metric][0], ci95_high=intervals[metric][1],
                pairs=len(frame), positives=int(y.sum()), clusters=frame.cluster_id.nunique(),
                precision_k=e['precision_k'], bootstrap_valid_draws=draws, scoring_seconds=round(score_seconds, 3)))
        frame['sequence_only_score'] = scores
        predictions.append(frame[['experimental_pair_id', 'vhh_id', 'cluster_id', 'split', 'Ag_label', 'label', 'sequence_only_score']])
        precision, recall, _ = precision_recall_curve(y, scores, drop_intermediate=True)
        ax.step(recall[::-1], precision[::-1], where='post', color='#147d92', linewidth=1.4,
                label=f'AP = {values["average_precision"]:.3f}')
        ax.axhline(values['prevalence'], color='#a65c16', linestyle='--',
                   label=f'Prevalence = {values["prevalence"]:.3f}')
        ax.set(title=f'{split.title()} ({len(frame):,} pairs)', xlabel='Recall', xlim=(0, 1), ylim=(0, 1))
        ax.legend(loc='upper right', fontsize=8)
        print(f'{split}: {values}', flush=True)
    axes[0].set_ylabel('Precision')
    fig.suptitle('AVIDa-hIL6 sequence-only baseline • held-out VHH clusters')
    fig.text(.5, .01, 'Original split prevalence; known IL-6 variants; no structural scores', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .035, 1, .94))
    fig.savefig(out / 'avida_sequence_baseline_pr.png', dpi=180)
    fig.savefig(out / 'avida_sequence_baseline_pr.pdf')
    plt.close(fig)
    pd.DataFrame(rows).to_csv(out / 'avida_sequence_baseline_metrics.csv', index=False)
    pd.concat(predictions).to_csv(ROOT / 'data/processed/avida_sequence_predictions.csv.gz', index=False,
                                  compression={'method':'gzip', 'mtime':0})
    metadata = dict(protocol_sha256=digest(cfg_path), train_pairs=len(train), train_prevalence=float(train.label.mean()),
        vocabulary=len(vectorizer.vocabulary_), antigen_categories=encoder.categories_[0].tolist(),
        sklearn_version=sklearn.__version__, numpy_version=np.__version__, model_iterations=int(model.n_iter_),
        fit_seconds=round(fit_seconds, 3), total_seconds=round(time.monotonic()-started, 3),
        warnings=[str(w.message) for w in caught], source_manifest_sha256=digest(out / 'avida_source_manifest.json'),
        model_sha256=digest(artifacts / 'model.npz'), predictions_sha256=digest(ROOT / 'data/processed/avida_sequence_predictions.csv.gz'),
        interpretation='Sequence-only baseline on released pair prevalence; not a structural-channel comparison or population-prevalence estimate',
        bootstrap='1000 component-resampling draws; fixed fitted model; percentile 95% intervals; variants remain together',
        tie_policy='Score descending then pair-ID ascending; same ordering in cluster bootstrap',
        test_policy='No hyperparameter selection or calibration on development/test; protocol fixed before fitting')
    (out / 'avida_baseline_run.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == '__main__':
    main()
