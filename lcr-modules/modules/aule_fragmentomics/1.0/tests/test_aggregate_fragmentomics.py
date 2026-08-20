#!/usr/bin/env python3

import importlib.util
import os
import tempfile
import unittest


HERE = os.path.dirname(os.path.abspath(__file__))
SOURCE = os.path.join(HERE, "..", "src", "aggregate_fragmentomics.py")
SPEC = importlib.util.spec_from_file_location("aggregate_fragmentomics", SOURCE)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class AggregateTests(unittest.TestCase):
    def test_aggregate_sorts_samples(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = []
            for sample in ("B", "A"):
                path = os.path.join(directory, sample + ".tsv")
                with open(path, "w") as handle:
                    handle.write("sample_id\tn_fragments\n{}\t1\n".format(sample))
                paths.append(path)
            rows, fields = MODULE.aggregate(paths)
            self.assertEqual(fields, ["sample_id", "n_fragments"])
            self.assertEqual([row["sample_id"] for row in rows], ["A", "B"])

    def test_duplicate_sample_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "duplicates.tsv")
            with open(path, "w") as handle:
                handle.write("sample_id\tn_fragments\nA\t1\nA\t2\n")
            with self.assertRaisesRegex(ValueError, "duplicate sample_id"):
                MODULE.aggregate([path])


    def test_regional_summary_is_joined_by_sample(self):
        with tempfile.TemporaryDirectory() as directory:
            global_path = os.path.join(directory, "global.tsv")
            regional_path = os.path.join(directory, "regional.tsv")
            with open(global_path, "w") as handle:
                handle.write("sample_id\tn_fragments\nA\t10\n")
            with open(regional_path, "w") as handle:
                handle.write("sample_id\tdelfi_status\nA\tpass\n")
            rows, fields = MODULE.aggregate([global_path], [regional_path])
            self.assertEqual(rows[0]["delfi_status"], "pass")
            self.assertIn("delfi_status", fields)


if __name__ == "__main__":
    unittest.main()
