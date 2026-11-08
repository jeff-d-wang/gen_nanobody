"""CPU-only checks for benchmark bookkeeping and output validation."""
import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import benchmark_rfantibody as bench


class BenchmarkChecks(unittest.TestCase):
    def test_prepared_target_mapping(self):
        rows = json.loads((bench.PROC / 'rfab_inputs/targets.json').read_text())
        self.assertEqual(len(rows), 21)
        self.assertEqual(rows[0]['hotspots'], ['T31', 'T38', 'T43', 'T47'])
        for row in rows:
            lines = Path(row['target']).read_text().splitlines()
            atoms = [s for s in lines if s.startswith('ATOM')]
            self.assertEqual({s[21] for s in atoms}, {'T'})
            ca = {int(s[22:26]) for s in atoms if s[12:16].strip() == 'CA'}
            self.assertEqual(len(ca), row['antigen_residues'])
            self.assertTrue({int(h[1:]) for h in row['hotspots']} <= ca)
            self.assertEqual(bench.digest(Path(row['target'])), row['target_sha256'])
        mse = next(r for r in rows if r['pdb_id'] == 'pdb_00008rwf')
        self.assertEqual(mse['mse_to_met_residues'], [85])
        self.assertEqual(mse['antigen_residues'], 165)

    def test_report_counts_failed_attempt_cost_and_partial_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'results').mkdir()
            run = root / 'run'
            run.mkdir()
            state = dict(pdb_id='test', arm='hotspot', status='failed', stages={},
                         attempts=[dict(stage='rfdiffusion', seconds=12, partial_outputs=[]),
                                   dict(stage='rfdiffusion', seconds=30, partial_outputs=[{}])])
            with patch.object(bench, 'ROOT', root):
                bench.report(run, [state])
            result = json.loads((run / 'summary.json').read_text())
            self.assertEqual(result['wall_seconds_all_attempts'], 42)
            self.assertEqual(result['runs'][0]['backbones'], 1)
            self.assertEqual(result['complete'], 0)

    def test_partial_stage_not_complete(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, 'expected 2'):
                bench.validate(Path(tmp), 'rfdiffusion')


if __name__ == '__main__':
    unittest.main()
