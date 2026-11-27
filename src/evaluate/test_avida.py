"""Leakage, label, weighted-bootstrap, and saved-metric regression tests."""
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src/data'))
from prepare_avida import components, split_components, validate_sources
import avida_baseline as baseline


class CurationTests(unittest.TestCase):
    def source(self):
        pairs = pd.DataFrame([dict(VHH_sequence='ACDEFGHIKL', Ag_label='WT', label=1,
                                   subject_species='alpaca', subject_name='test', subject_sex='male')])
        return pairs, pd.DataFrame([dict(Ag_label='WT', Ag_sequence='ACDEFGHIKLMNPQ')])

    def test_transitive_components_and_indel_boundary(self):
        seqs = ['AAAAAAAAAA', 'CAAAAAAAAA', 'CCAAAAAAAA', 'GGGGGGGG', 'GGGGGGGGG', 'GGGGGGGGGG']
        groups, stats = components(seqs, block_size=2, workers=1, progress=False)
        self.assertEqual(len(set(groups[:3])), 1)  # endpoints differ by 20%, but bridge stays together
        self.assertNotEqual(groups[3], groups[4])  # 1/9 exceeds 10%
        self.assertEqual(groups[4], groups[5])     # 1/10 is accepted
        self.assertEqual(stats['pairwise_comparisons'], 15)

    def test_label_conflicts_rejected(self):
        pairs, antigens = self.source()
        conflict = pairs.copy()
        conflict['label'] = 0
        with self.assertRaisesRegex(ValueError, 'Conflicting'):
            validate_sources(pd.concat([pairs, conflict], ignore_index=True), antigens)

    def test_unmatched_antigen_and_unknown_amino_acid_rejected(self):
        pairs, antigens = self.source()
        antigens['Ag_label'] = 'OTHER'
        with self.assertRaisesRegex(ValueError, 'Unmatched'):
            validate_sources(pairs, antigens)
        pairs, antigens = self.source()
        pairs['VHH_sequence'] = 'ACDX'
        with self.assertRaisesRegex(ValueError, 'Noncanonical'):
            validate_sources(pairs, antigens)

    def test_exact_duplicates_counted_and_collapsed(self):
        pairs, antigens = self.source()
        joined, removed = validate_sources(pd.concat([pairs, pairs], ignore_index=True), antigens)
        self.assertEqual((len(joined), removed), (1, 1))

    def test_split_independent_of_input_order(self):
        frame = pd.DataFrame(dict(vhh_id=list('abcdef'), cluster_id=['A', 'A', 'B', 'C', 'D', 'E']))
        counts = dict(zip('abcdef', [10, 2, 8, 4, 1, 3]))
        first = frame.assign(split=split_components(frame, counts)).set_index('vhh_id').split
        shuffled = frame.sample(frac=1, random_state=11)
        second = shuffled.assign(split=split_components(shuffled, counts)).set_index('vhh_id').split
        self.assertEqual(first.to_dict(), second.to_dict())
        self.assertEqual(first.a, first.b)


class MetricTests(unittest.TestCase):
    def test_weighted_bootstrap_equals_explicit_cluster_replication(self):
        y = np.array([1, 0, 1, 0, 0])
        scores = np.array([.8, .8, .7, .3, .1])
        weights = np.array([2, 2, 0, 1, 3])
        weighted = baseline.metric_values(y, scores, 4, weights)
        expanded = baseline.metric_values(np.repeat(y, weights), np.repeat(scores, weights), 4)
        for key in weighted:
            self.assertAlmostEqual(weighted[key], expanded[key])

    def test_bootstrap_keeps_variants_together(self):
        captured = []
        def metric(y, score, k, weights):
            captured.append(weights)
            return dict(prevalence=float(np.average(y, weights=weights)))
        with patch.object(baseline, 'metric_values', metric):
            baseline.bootstrap(np.array([0, 1, 0, 1]), np.arange(4), np.array(['a', 'a', 'b', 'b']), 2, 10, 37)
        for weights in captured:
            self.assertEqual(weights[0], weights[1])
            self.assertEqual(weights[2], weights[3])

    def test_saved_test_metrics_independently_recomputed(self):
        predictions = pd.read_csv(ROOT / 'data/processed/avida_sequence_predictions.csv.gz')
        frame = predictions[predictions.split == 'test'].sort_values('experimental_pair_id')
        y, scores = frame.label.to_numpy(), frame.sequence_only_score.to_numpy()
        order = np.argsort(-scores, kind='stable')
        cumulative = np.cumsum(y[order])
        ends = np.r_[np.flatnonzero(np.diff(scores[order]) != 0), len(y)-1]
        ap = np.sum(np.diff(np.r_[0, cumulative[ends]]) / y.sum() * cumulative[ends] / (ends+1))
        positives, negatives = y.sum(), len(y)-y.sum()
        auc = (rankdata(scores)[y == 1].sum() - positives*(positives+1)/2) / (positives*negatives)
        expected = dict(average_precision=ap, auroc=auc, precision_at_k=y[order[:100]].mean(), prevalence=y.mean())
        saved = pd.read_csv(ROOT / 'results/evaluation/avida_sequence_baseline_metrics.csv')
        for metric, value in expected.items():
            actual = saved.loc[(saved.split == 'test') & (saved.metric == metric), 'value'].item()
            self.assertAlmostEqual(value, actual, places=10)

    def test_saved_splits_and_prediction_pair_identity(self):
        pairs = pd.read_csv(ROOT / 'data/processed/avida_pairs.csv.gz')
        self.assertTrue(pairs.experimental_pair_id.is_unique)
        self.assertEqual(pairs.groupby('cluster_id').split.nunique().max(), 1)
        self.assertEqual(pairs.groupby('vhh_id').split.nunique().max(), 1)
        predictions = pd.read_csv(ROOT / 'data/processed/avida_sequence_predictions.csv.gz')
        joined = predictions.merge(pairs[['experimental_pair_id', 'label', 'split']], on='experimental_pair_id', validate='one_to_one')
        self.assertEqual(len(joined), len(predictions))
        self.assertTrue((joined.label_x == joined.label_y).all())
        self.assertTrue((joined.split_x == joined.split_y).all())
        self.assertEqual(set(predictions.experimental_pair_id), set(pairs.loc[pairs.split != 'train', 'experimental_pair_id']))


if __name__ == '__main__':
    unittest.main()
