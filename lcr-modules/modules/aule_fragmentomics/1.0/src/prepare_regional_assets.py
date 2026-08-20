#!/usr/bin/env python3

"""Prepare deterministic DELFI windows and Griffin TF-site catalogues."""

import argparse
import gzip
import os
import sys
from collections import defaultdict


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fai(path):
    lengths = {}
    order = []
    with open(path) as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError("invalid FAI line {} in {}".format(line_number, path))
            chrom, length = fields[0], int(fields[1])
            if chrom in lengths:
                raise ValueError("duplicate contig {} in {}".format(chrom, path))
            lengths[chrom] = length
            order.append(chrom)
    return lengths, order


def read_bed_intervals(path, known_contigs=None):
    intervals = defaultdict(list)
    if not path:
        return intervals
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError("BED line {} has fewer than three columns".format(line_number))
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError("invalid BED interval on line {}".format(line_number))
            if known_contigs is not None and chrom not in known_contigs:
                continue
            intervals[chrom].append((start, end))
    for chrom in intervals:
        intervals[chrom] = merge_intervals(intervals[chrom])
    return intervals


def merge_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [tuple(interval) for interval in merged]


def overlap_length(intervals, start, end):
    total = 0
    for left, right in intervals:
        if right <= start:
            continue
        if left >= end:
            break
        total += max(0, min(end, right) - max(start, left))
    return total


def parse_chromosomes(value, fai_order):
    if not value or value.lower() == "autosomes":
        expected = ["chr{}".format(number) for number in range(1, 23)]
        if not any(chrom in fai_order for chrom in expected):
            expected = [str(number) for number in range(1, 23)]
        return [chrom for chrom in expected if chrom in fai_order]
    requested = [item.strip() for item in value.split(",") if item.strip()]
    missing = [chrom for chrom in requested if chrom not in fai_order]
    if missing:
        raise ValueError("contigs absent from FAI: {}".format(", ".join(missing)))
    return requested


def gc_fraction(sequence):
    sequence = sequence.upper()
    acgt = sum(sequence.count(base) for base in "ACGT")
    if acgt == 0:
        return None
    return float(sequence.count("G") + sequence.count("C")) / acgt


def prepare_delfi(options):
    try:
        import pysam
    except ImportError as error:
        raise RuntimeError("pysam is required to generate DELFI bins") from error

    lengths, order = read_fai(options.fai)
    chromosomes = parse_chromosomes(options.chromosomes, order)
    if not chromosomes:
        raise ValueError("no selected chromosomes are present in the FAI")
    blacklist = read_bed_intervals(options.blacklist, lengths)
    os.makedirs(os.path.dirname(os.path.abspath(options.output)), exist_ok=True)
    written = 0
    skipped = 0
    fasta = pysam.FastaFile(options.fasta)
    try:
        with open(options.output, "w") as output:
            for chrom in chromosomes:
                for start in range(0, lengths[chrom], options.bin_size):
                    end = min(start + options.bin_size, lengths[chrom])
                    blocked = overlap_length(blacklist.get(chrom, []), start, end)
                    blocked_fraction = float(blocked) / (end - start)
                    if blocked_fraction > options.max_blacklist_fraction:
                        skipped += 1
                        continue
                    gc = gc_fraction(fasta.fetch(chrom, start, end))
                    if gc is None:
                        skipped += 1
                        continue
                    bin_id = "{}:{:09d}-{:09d}".format(chrom, start, end)
                    output.write("{}\t{}\t{}\t{}\t{:.8f}\t{:.8f}\n".format(
                        chrom, start, end, bin_id, gc, blocked_fraction))
                    written += 1
    finally:
        fasta.close()
    if written == 0:
        raise ValueError("DELFI preparation produced no bins")
    print("DELFI bins written: {}; filtered: {}".format(written, skipped), file=sys.stderr)


def read_allowlist(path):
    if not path:
        return None
    result = set()
    with open_text(path) as handle:
        for line in handle:
            value = line.strip().split("\t")[0]
            if value and not value.startswith("#"):
                result.add(value)
    return result


def prepare_griffin(options):
    lengths, order = read_fai(options.fai)
    rank = {chrom: index for index, chrom in enumerate(order)}
    blacklist = read_bed_intervals(options.blacklist, lengths)
    allowlist = read_allowlist(options.tf_allowlist)
    rows = {}
    filtered = defaultdict(int)
    with open_text(options.input) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError("Griffin BED line {} requires at least BED4".format(line_number))
            chrom, start, end, tf = fields[0], int(fields[1]), int(fields[2]), fields[3]
            if chrom not in lengths:
                filtered["unknown_contig"] += 1
                continue
            if start < 0 or end <= start or end > lengths[chrom]:
                raise ValueError("invalid Griffin interval on line {}".format(line_number))
            if allowlist is not None and tf not in allowlist:
                filtered["not_allowlisted"] += 1
                continue
            centre = (start + end) // 2
            if centre - options.window < 0 or centre + options.window > lengths[chrom]:
                filtered["contig_edge"] += 1
                continue
            if overlap_length(blacklist.get(chrom, []), start, end) > 0:
                filtered["blacklist"] += 1
                continue
            key = (chrom, start, end, tf)
            site_id = fields[4] if len(fields) >= 5 and fields[4] else \
                "{}:{}-{}:{}".format(chrom, start, end, tf)
            if key in rows:
                filtered["duplicate"] += 1
                continue
            rows[key] = site_id
    sorted_rows = sorted(rows, key=lambda row: (rank[row[0]], row[1], row[2], row[3]))
    os.makedirs(os.path.dirname(os.path.abspath(options.output)), exist_ok=True)
    with open(options.output, "w") as output:
        for chrom, start, end, tf in sorted_rows:
            output.write("{}\t{}\t{}\t{}\t{}\n".format(
                chrom, start, end, tf, rows[(chrom, start, end, tf)]))
    if not sorted_rows:
        raise ValueError("Griffin preparation produced no sites")
    details = ", ".join("{}={}".format(key, filtered[key]) for key in sorted(filtered))
    print("Griffin sites written: {}{}".format(
        len(sorted_rows), "; filtered: " + details if details else ""), file=sys.stderr)


def parser():
    result = argparse.ArgumentParser()
    subparsers = result.add_subparsers(dest="command", required=True)
    delfi = subparsers.add_parser("delfi-bins")
    delfi.add_argument("--fasta", required=True)
    delfi.add_argument("--fai", required=True)
    delfi.add_argument("--blacklist", required=True)
    delfi.add_argument("--output", required=True)
    delfi.add_argument("--bin-size", type=int, default=5000000)
    delfi.add_argument("--chromosomes", default="autosomes")
    delfi.add_argument("--max-blacklist-fraction", type=float, default=0.10)
    delfi.set_defaults(function=prepare_delfi)
    griffin = subparsers.add_parser("griffin-sites")
    griffin.add_argument("--input", required=True)
    griffin.add_argument("--fai", required=True)
    griffin.add_argument("--blacklist", required=True)
    griffin.add_argument("--output", required=True)
    griffin.add_argument("--window", type=int, default=1000)
    griffin.add_argument("--tf-allowlist", default="")
    griffin.set_defaults(function=prepare_griffin)
    return result


def main(argv=None):
    options = parser().parse_args(argv)
    if getattr(options, "bin_size", 1) <= 0:
        raise ValueError("bin-size must be positive")
    if not 0 <= getattr(options, "max_blacklist_fraction", 0) <= 1:
        raise ValueError("max-blacklist-fraction must be in [0, 1]")
    if getattr(options, "window", 0) < 0:
        raise ValueError("window must be non-negative")
    options.function(options)


if __name__ == "__main__":
    main()
