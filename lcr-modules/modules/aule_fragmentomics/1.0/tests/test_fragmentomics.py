#!/usr/bin/env python3

import importlib.util
import math
import os
import tempfile
import unittest
from collections import Counter


HERE = os.path.dirname(os.path.abspath(__file__))
SOURCE = os.path.join(HERE, "..", "src", "collect_fragmentomics.py")
SPEC = importlib.util.spec_from_file_location("collect_fragmentomics", SOURCE)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class Options:
    min_length = 50
    max_length = 250
    short_min = 100
    short_max = 150
    long_min = 151
    long_max = 220
    periodicity_max_length = 150
    periodicity_min_bp = 8
    periodicity_max_bp = 12


class FakeRead:
    is_proper_pair = True
    mapping_quality = 60
    flag = 0
    query_sequence = "A"

    def __init__(self, template_length, reference_start, reference_end, quality=30):
        self.template_length = template_length
        self.reference_start = reference_start
        self.reference_end = reference_end
        self.query_qualities = [quality]

    def get_aligned_pairs(self, matches_only=False):
        return [(0, 100)]


class FragmentomicsTests(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(MODULE.reverse_complement("ACGT"), "ACGT")
        self.assertEqual(MODULE.reverse_complement("AAAC"), "GTTT")

    def test_weighted_median_and_size_metrics(self):
        histogram = Counter({100: 1, 110: 2, 160: 1})
        self.assertEqual(MODULE.weighted_median(histogram), 110)
        metrics = MODULE.size_metrics(histogram, Options())
        self.assertEqual(metrics["n_fragments"], 4)
        self.assertEqual(metrics["modal_fragment_length"], 110)
        self.assertEqual(metrics["p_short"], 0.75)
        self.assertEqual(metrics["short_long_ratio"], 3.0)

    def test_motif_diversity_bounds(self):
        self.assertEqual(MODULE.normalised_entropy(Counter({"AAAA": 5}), 256), 0.0)
        uniform = Counter({motif: 1 for motif in MODULE.MOTIFS})
        self.assertAlmostEqual(MODULE.normalised_entropy(uniform, 256), 1.0)

    def test_periodicity_detects_ten_bp_pattern(self):
        histogram = Counter({length: 100 for length in range(50, 151, 10)})
        amplitude, period = MODULE.periodicity_score(histogram, 50, 150, 8, 12)
        self.assertEqual(period, 10)
        self.assertGreater(amplitude, 0.5)

    def test_ks_identical_and_separated(self):
        distance, p_value = MODULE.ks_2sample([100, 110, 120], [100, 110, 120])
        self.assertEqual(distance, 0)
        self.assertEqual(p_value, 1)
        distance, p_value = MODULE.ks_2sample([80] * 20, [200] * 20)
        self.assertEqual(distance, 1)
        self.assertLess(p_value, 1e-6)

    def test_incomplete_weight_table_is_not_partially_applied(self):
        result = MODULE.apply_weights(Counter({100: 2}), Counter({110: 3}), {100: 2.0})
        self.assertEqual(result["weighting_status"], "incomplete_weight_table")
        self.assertIsNone(result["weighted_variant_vaf"])

    def test_weight_table_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "weights.tsv")
            with open(path, "w") as handle:
                handle.write("length\tlikelihood_ratio\n100\t2\n110\t0.5\n")
            self.assertEqual(MODULE.load_size_weights(path), {100: 2.0, 110: 0.5})

    def test_allele_filter_accepts_negative_tlen_mate(self):
        read = FakeRead(-180, 180, 230)
        self.assertTrue(MODULE.allele_read_passes(read, 30, 3852, 50, 250))
        self.assertEqual(MODULE.fragment_bounds(read), (50, 230))

    def test_variant_base_requires_minimum_base_quality(self):
        low_quality = FakeRead(180, 50, 101, quality=29)
        passing = FakeRead(180, 50, 101, quality=30)
        self.assertIsNone(MODULE.base_at_reference_position(low_quality, 100, 30))
        self.assertEqual(MODULE.base_at_reference_position(passing, 100, 30), "A")


if __name__ == "__main__":
    unittest.main()
