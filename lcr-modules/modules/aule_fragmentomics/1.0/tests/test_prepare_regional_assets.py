import importlib.util
import os
import sys
import tempfile
import types
import unittest


HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(HERE), "src", "prepare_regional_assets.py")
SPEC = importlib.util.spec_from_file_location("prepare_regional_assets", SCRIPT)
ASSETS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ASSETS)


class FakeFasta:
    sequences = {}

    def __init__(self, path):
        self.path = path

    def fetch(self, chrom, start, end):
        return self.sequences[chrom][start:end]

    def close(self):
        pass


class PrepareRegionalAssetsTests(unittest.TestCase):
    def test_merge_and_overlap_do_not_double_count_blacklist(self):
        merged = ASSETS.merge_intervals([(8, 15), (2, 10), (20, 22)])
        self.assertEqual(merged, [(2, 15), (20, 22)])
        self.assertEqual(ASSETS.overlap_length(merged, 5, 21), 11)

    def test_delfi_bins_include_gc_and_blacklist_fraction(self):
        with tempfile.TemporaryDirectory() as directory:
            fai = os.path.join(directory, "reference.fa.fai")
            blacklist = os.path.join(directory, "blacklist.bed")
            output = os.path.join(directory, "delfi.bed")
            with open(fai, "w") as handle:
                handle.write("chr1\t20\t0\t0\t0\n")
            with open(blacklist, "w") as handle:
                handle.write("chr1\t2\t4\n")
                handle.write("chr1\t10\t18\n")
            FakeFasta.sequences = {"chr1": "ACGT" * 5}
            previous = sys.modules.get("pysam")
            sys.modules["pysam"] = types.SimpleNamespace(FastaFile=FakeFasta)
            try:
                options = types.SimpleNamespace(
                    fasta="unused.fa", fai=fai, blacklist=blacklist,
                    output=output, bin_size=10, chromosomes="chr1",
                    max_blacklist_fraction=0.5)
                ASSETS.prepare_delfi(options)
            finally:
                if previous is None:
                    del sys.modules["pysam"]
                else:
                    sys.modules["pysam"] = previous
            with open(output) as handle:
                rows = [line.rstrip().split("\t") for line in handle]
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0][:4], ["chr1", "0", "10", "chr1:000000000-000000010"])
            self.assertEqual(rows[0][4:], ["0.50000000", "0.20000000"])

    def test_griffin_sites_are_filtered_deduplicated_and_sorted(self):
        with tempfile.TemporaryDirectory() as directory:
            fai = os.path.join(directory, "reference.fa.fai")
            blacklist = os.path.join(directory, "blacklist.bed")
            raw = os.path.join(directory, "raw.bed")
            allowlist = os.path.join(directory, "tf.txt")
            output = os.path.join(directory, "sites.bed")
            with open(fai, "w") as handle:
                handle.write("chr1\t10000\t0\t0\t0\nchr2\t10000\t0\t0\t0\n")
            with open(blacklist, "w") as handle:
                handle.write("chr1\t4900\t5100\n")
            with open(allowlist, "w") as handle:
                handle.write("TF1\n")
            with open(raw, "w") as handle:
                handle.write("chr2\t3000\t3010\tTF1\tsite2\n")
                handle.write("chr1\t3000\t3010\tTF1\tsite1\n")
                handle.write("chr1\t3000\t3010\tTF1\tduplicate\n")
                handle.write("chr1\t5000\t5010\tTF1\tblacklisted\n")
                handle.write("chr1\t100\t110\tTF1\tedge\n")
                handle.write("chr1\t4000\t4010\tTF2\tnot-allowed\n")
            options = types.SimpleNamespace(
                input=raw, fai=fai, blacklist=blacklist, output=output,
                window=1000, tf_allowlist=allowlist)
            ASSETS.prepare_griffin(options)
            with open(output) as handle:
                rows = [line.rstrip() for line in handle]
            self.assertEqual(rows, [
                "chr1\t3000\t3010\tTF1\tsite1",
                "chr2\t3000\t3010\tTF1\tsite2",
            ])


if __name__ == "__main__":
    unittest.main()
