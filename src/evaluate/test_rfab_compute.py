"""Regression checks for retry accounting and the frozen archive's independent FASTA."""
import csv
import json
import unittest

from rfab_compute import ROOT, STAGES, attempt_rows


class ComputeTests(unittest.TestCase):
    def state(self):
        attempts = []
        for stage, number, seconds in [('rfdiffusion', 1, 17), ('rfdiffusion', 2, 23),
                                       ('proteinmpnn', 1, 3), ('rf2', 1, 7)]:
            attempts.append(dict(stage=stage, log=f'/old/machine/{stage}_attempt_{number}/output.log',
                seconds=seconds, started_utc='2026-11-07T00:00:00+00:00',
                finished_utc='2026-11-07T00:01:00+00:00', returncode=0))
        accepted = {s: f'{s}_attempt_{2 if s == "rfdiffusion" else 1}' for s in STAGES}
        return dict(pdb_id='test', arm='unguided', attempts=attempts), accepted

    def test_retry_cost_included_without_duplicate_yield(self):
        state, accepted = self.state()
        rows = attempt_rows(state, accepted)
        self.assertEqual(sum(r['wall_seconds'] for r in rows), 50)
        self.assertEqual(sum(r['wall_seconds'] for r in rows if not r['accepted']), 17)
        self.assertEqual(sum(r['audited_candidates_added'] for r in rows), 2)

    def test_missing_or_nonfinite_cost_rejected(self):
        for value in (float('nan'), float('inf'), -1):
            state, accepted = self.state()
            state['attempts'][0]['seconds'] = value
            with self.assertRaises(ValueError):
                attempt_rows(state, accepted)
        state, accepted = self.state()
        del state['attempts'][0]['seconds']
        with self.assertRaises(KeyError):
            attempt_rows(state, accepted)

    def test_missing_accepted_attempt_rejected(self):
        state, accepted = self.state()
        accepted['rf2'] = 'rf2_attempt_9'
        with self.assertRaises(ValueError):
            attempt_rows(state, accepted)

    def test_archived_fasta_independently_matches_structure_sequences(self):
        run = ROOT / 'runs/rfab_benchmark_1ce3344e4822'
        aa = dict(pair.split(':') for pair in
            'ALA:A ARG:R ASN:N ASP:D CYS:C GLN:Q GLU:E GLY:G HIS:H ILE:I LEU:L LYS:K MET:M PHE:F PRO:P SER:S THR:T TRP:W TYR:Y VAL:V'.split())
        fasta = {}
        for line in (run / 'designed_sequences.fasta').read_text().splitlines():
            if line.startswith('>'):
                key = line[1:]
                self.assertNotIn(key, fasta)
                fasta[key] = ''
            else:
                fasta[key] += line.strip()
        with (ROOT / 'results/evaluation/rfab_compute_candidates.csv').open() as stream:
            candidates = list(csv.DictReader(stream))
        observed = {}
        for row in candidates:
            key = f'{row["pdb_id"]}_{row["arm"]}_nb_{row["backbone_index"]}_dldesign_0'
            observed[key] = ''.join(aa[line[17:20]] for line in (ROOT / row['designed']).read_text().splitlines()
                                   if line.startswith('ATOM') and line[21] == 'H' and line[12:16].strip() == 'CA')
        self.assertEqual(observed, fasta)
        self.assertEqual(len(fasta), 84)
        self.assertEqual(len(set(fasta.values())), 83)


if __name__ == '__main__':
    unittest.main()
