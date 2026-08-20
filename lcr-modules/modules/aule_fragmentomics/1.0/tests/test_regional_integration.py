#!/usr/bin/env python3

import csv
import importlib.util
import os
import tempfile
import unittest
from types import SimpleNamespace

try:
    import pysam
except ImportError:  # The module conda environment supplies pysam.
    pysam = None


HERE = os.path.dirname(os.path.abspath(__file__))
SOURCE = os.path.join(HERE, "..", "src", "collect_regional_fragmentomics.py")
SPEC = importlib.util.spec_from_file_location("collect_regional_integration", SOURCE)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


@unittest.skipIf(pysam is None, "pysam is available in the module environment")
class RegionalIntegrationTests(unittest.TestCase):
    def test_synthetic_bam_produces_delfi_and_griffin_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta = os.path.join(directory, "reference.fa")
            with open(fasta, "w") as handle:
                handle.write(">chr1\n" + ("AAAACCCCGGGGTTTT" * 250) + "\n")
            pysam.faidx(fasta)

            unsorted = os.path.join(directory, "reads.unsorted.bam")
            bam = os.path.join(directory, "reads.bam")
            header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 4000}]}
            with pysam.AlignmentFile(unsorted, "wb", header=header) as output:
                read_number = 0
                for centre in (500, 1500, 2500):
                    for length in (120, 130, 160, 180):
                        start = centre - length // 2
                        left = pysam.AlignedSegment()
                        left.query_name = "fragment{}".format(read_number)
                        left.query_sequence = "A" * 50
                        left.flag = 99
                        left.reference_id = 0
                        left.reference_start = start
                        left.mapping_quality = 60
                        left.cigar = ((0, 50),)
                        left.next_reference_id = 0
                        left.next_reference_start = start + length - 50
                        left.template_length = length
                        left.query_qualities = pysam.qualitystring_to_array("I" * 50)
                        right = pysam.AlignedSegment()
                        right.query_name = left.query_name
                        right.query_sequence = "T" * 50
                        right.flag = 147
                        right.reference_id = 0
                        right.reference_start = start + length - 50
                        right.mapping_quality = 60
                        right.cigar = ((0, 50),)
                        right.next_reference_id = 0
                        right.next_reference_start = start
                        right.template_length = -length
                        right.query_qualities = pysam.qualitystring_to_array("I" * 50)
                        output.write(left)
                        output.write(right)
                        read_number += 1
            pysam.sort("-o", bam, unsorted)
            pysam.index(bam)

            blacklist = os.path.join(directory, "blacklist.bed")
            open(blacklist, "w").close()
            bins = os.path.join(directory, "bins.bed")
            with open(bins, "w") as handle:
                handle.write("chr1\t0\t1000\tbin1\nchr1\t1000\t2000\tbin2\n"
                             "chr1\t2000\t3000\tbin3\n")
            sites = os.path.join(directory, "sites.bed")
            with open(sites, "w") as handle:
                handle.write("chr1\t499\t501\tPAX5\ts1\n"
                             "chr1\t1499\t1501\tPAX5\ts2\n"
                             "chr1\t2499\t2501\tPAX5\ts3\n")

            options = SimpleNamespace(
                sample_id="S1", assay_partition="off_target", bam=bam, fasta=fasta, blacklist=blacklist,
                enable_delfi=True, delfi_bins=bins, delfi_reference="",
                enable_griffin=True, griffin_sites=sites, griffin_reference="",
                min_mapq=30, exclude_flags=3852, min_length=50, max_length=250,
                short_min=100, short_max=150, long_min=151, long_max=220,
                delfi_pseudocount=0.5, delfi_min_fragments=2,
                delfi_min_bins_for_lowess=3, delfi_gc_span=1.0,
                griffin_window_bp=100, griffin_bin_size_bp=10,
                griffin_central_bp=30, griffin_amplitude_bp=100,
                griffin_min_bins_for_lowess=3, griffin_min_sites=3,
                griffin_gc_span=1.0,
                delfi_output=os.path.join(directory, "delfi.tsv"),
                griffin_profile_output=os.path.join(directory, "profile.tsv"),
                griffin_features_output=os.path.join(directory, "features.tsv"),
                summary_output=os.path.join(directory, "summary.tsv"))
            MODULE.run(options)

            with open(options.delfi_output, newline="") as handle:
                delfi_rows = list(csv.DictReader(handle, delimiter="\t"))
            with open(options.griffin_profile_output, newline="") as handle:
                profile_rows = list(csv.DictReader(handle, delimiter="\t"))
            with open(options.griffin_features_output, newline="") as handle:
                feature_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(delfi_rows), 3)
            self.assertTrue(all(row["gc_corrected_log2_ratio"] != "NA"
                                for row in delfi_rows))
            self.assertEqual(len(profile_rows), 20)
            self.assertEqual({row["feature"] for row in feature_rows},
                             {"central_coverage", "mean_coverage", "amplitude"})
            self.assertTrue(all(row["status"] == "pass" for row in feature_rows))


if __name__ == "__main__":
    unittest.main()
