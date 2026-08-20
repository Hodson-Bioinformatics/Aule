#!/usr/bin/env python3

import csv
import importlib.util
import os
import tempfile
import unittest
from collections import Counter
from types import SimpleNamespace


HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, filename):
    path = os.path.join(HERE, "..", "src", filename)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


REGIONAL = load("collect_regional_fragmentomics", "collect_regional_fragmentomics.py")
REFERENCE = load("build_regional_reference", "build_regional_reference.py")


class RegionalTests(unittest.TestCase):
    def test_lowess_recovers_linear_gc_effect(self):
        x = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65]
        y = [1 + 2 * value for value in x]
        fitted = REGIONAL.lowess(x, y, span=0.75)
        for observed, expected in zip(fitted, y):
            self.assertAlmostEqual(observed, expected, places=8)

    def test_delfi_bins_use_fragment_midpoint(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "bins.bed")
            with open(path, "w") as handle:
                handle.write("chr1\t0\t100\tbin1\nchr1\t100\t200\tbin2\n")
            bins = REGIONAL.DelfiBins(path)
            self.assertEqual(bins.midpoint_bin("chr1", 99)["bin_id"], "bin1")
            self.assertEqual(bins.midpoint_bin("chr1", 100)["bin_id"], "bin2")
            self.assertIsNone(bins.midpoint_bin("chr1", 200))

    def test_overlapping_delfi_bins_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "bins.bed")
            with open(path, "w") as handle:
                handle.write("chr1\t0\t101\ta\nchr1\t100\t200\tb\n")
            with self.assertRaises(ValueError):
                REGIONAL.DelfiBins(path)

    def test_griffin_half_open_profile_bins(self):
        self.assertEqual(REGIONAL.griffin_bin_bounds(-1000, -990, 1000, 10), (0, 1))
        self.assertEqual(REGIONAL.griffin_bin_bounds(-5, 5, 1000, 10), (99, 101))
        self.assertEqual(REGIONAL.griffin_bin_bounds(990, 1000, 1000, 10), (199, 200))

    def test_reference_builder_uses_independent_samples(self):
        values = {("bin1",): {"N1": 1.0, "N2": 2.0, "N3": 3.0}}
        rows = REFERENCE.reference_rows(values, ["bin_id"], 3)
        self.assertEqual(rows[0]["bin_id"], "bin1")
        self.assertEqual(rows[0]["reference_mean"], 2.0)
        self.assertEqual(rows[0]["reference_sd"], 1.0)
        self.assertEqual(rows[0]["n_reference_samples"], 3)

    def test_disabled_mode_does_not_open_bam(self):
        with tempfile.TemporaryDirectory() as directory:
            options = SimpleNamespace(
                enable_delfi=False, enable_griffin=False, sample_id="S1", assay_partition="off_target",
                bam="missing.bam", fasta="missing.fa", blacklist="missing.bed",
                delfi_output=os.path.join(directory, "delfi.tsv"),
                griffin_profile_output=os.path.join(directory, "profile.tsv"),
                griffin_features_output=os.path.join(directory, "features.tsv"),
                summary_output=os.path.join(directory, "summary.tsv"))
            REGIONAL.run(options)
            with open(options.summary_output, newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["delfi_status"], "not_requested")
            self.assertEqual(row["griffin_status"], "not_requested")


if __name__ == "__main__":
    unittest.main()
