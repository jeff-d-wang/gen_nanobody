"""Archive checks for the exact experimental-pair RF2 pilot (RFantibody env)."""
import csv
import json
import unittest

import numpy as np
from scipy.spatial import cKDTree

from avida_pilot import ROOT, target_template, inspect_prediction


class PilotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with (ROOT / 'results/evaluation/avida_pilot_pairs.csv').open() as stream:
            cls.pairs = {r['experimental_pair_id']: r for r in csv.DictReader(stream)}
        cls.seq = next(iter(cls.pairs.values()))['Ag_sequence']
        cls.template, cls.supported, cls.offset = target_template(ROOT / 'data/raw/avida_pilot/1ALU.pdb', cls.seq)

    def test_target_mapping_masks_unobserved_residues(self):
        self.assertEqual((len(self.seq), len(self.supported), self.offset), (218, 157, 27))
        xyz, mask = np.array(self.template['xyz']), np.array(self.template['mask'])
        self.assertTrue(np.isfinite(xyz[mask]).all())
        self.assertTrue(np.isnan(xyz[~mask]).all())
        self.assertFalse(mask[:27].any())
        self.assertFalse(mask[-6:].any())

    def test_template_sequence_mismatch_is_rejected(self):
        index = self.supported[0]
        bad = self.seq[:index] + ('A' if self.seq[index] != 'A' else 'G') + self.seq[index+1:]
        with self.assertRaisesRegex(ValueError, 'sequence mismatch'):
            target_template(ROOT / 'data/raw/avida_pilot/1ALU.pdb', bad)

    def test_all_saved_predictions_and_independent_geometry(self):
        with (ROOT / 'results/evaluation/avida_feasibility.csv').open() as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 4)
        self.assertEqual({r['experimental_pair_id'] for r in rows}, set(self.pairs))
        for row in rows:
            self.assertEqual(row['status'], 'passed')
            pair = self.pairs[row['experimental_pair_id']]
            path = ROOT / row['structure']
            values = inspect_prediction(path, pair['VHH_sequence'], pair['Ag_sequence'], self.supported)
            atoms = {'H': [], 'T': []}
            keys = {'H': [], 'T': []}
            for line in path.read_text().splitlines():
                if line.startswith('ATOM') and line[76:78].strip() != 'H' and not line[12:16].strip().startswith('H'):
                    atoms[line[21]].append([float(line[i:i+8]) for i in (30, 38, 46)])
                    keys[line[21]].append(line[22:27])
            tree_h, tree_t = cKDTree(atoms['H']), cKDTree(atoms['T'])
            neighbours = tree_h.query_ball_tree(tree_t, r=5.)
            contact_keys = {keys['T'][j] for js in neighbours for j in js}
            self.assertEqual(len(contact_keys), values['target_contact_residues_5A'])
            severe = sum(np.linalg.norm(np.array(atoms['H'][i])-atoms['T'][j]) < 2.
                         for i, js in enumerate(tree_h.query_ball_tree(tree_t, r=2.)) for j in js)
            self.assertEqual(severe, values['interchain_atom_pairs_below_2A'])
            for key, value in values.items():
                self.assertAlmostEqual(value, float(row[key]))
            with self.assertRaisesRegex(ValueError, 'changed chain H'):
                wrong = ('A' if pair['VHH_sequence'][0] != 'A' else 'G') + pair['VHH_sequence'][1:]
                inspect_prediction(path, wrong, pair['Ag_sequence'], self.supported)


if __name__ == '__main__':
    unittest.main()
