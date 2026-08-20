#!/usr/bin/env python3

"""DELFI-style and Griffin-style regional cfDNA fragmentomics.

The collector deliberately separates directly measured profiles from application
against an independently built healthy reference.  It does not emit a cancer
probability or train a classifier on the analysed cohort.
"""

import argparse
import bisect
import csv
import gzip
import math
import os
from collections import Counter, defaultdict


def parse_bool(value):
    value = str(value).strip().lower()
    if value in {"true", "1", "yes"}:
        return True
    if value in {"false", "0", "no"}:
        return False
    raise argparse.ArgumentTypeError("expected true or false")


def finite(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def display(value):
    if value is None or (not isinstance(value, (str, bool)) and not finite(value)):
        return "NA"
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    return value


def write_tsv(path, rows, fields):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields,
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: display(row.get(field)) for field in fields})


def bed_rows(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError("{} line {} has fewer than three columns".format(
                    path, line_number))
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError("{} line {} has invalid coordinates".format(
                    path, line_number))
            yield fields[0], start, end, fields[3:]


class Blacklist:
    def __init__(self, path):
        intervals = defaultdict(list)
        for chrom, start, end, _ in bed_rows(path):
            intervals[chrom].append((start, end))
        self.intervals = {}
        self.starts = {}
        for chrom, values in intervals.items():
            merged = []
            for start, end in sorted(values):
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    merged.append((start, end))
            self.intervals[chrom] = merged
            self.starts[chrom] = [start for start, _ in merged]

    def overlaps(self, chrom, start, end):
        values = self.intervals.get(chrom, [])
        if not values:
            return False
        index = bisect.bisect_left(self.starts[chrom], end) - 1
        return index >= 0 and values[index][1] > start


class DelfiBins:
    def __init__(self, path):
        self.rows = []
        by_chrom = defaultdict(list)
        seen = set()
        for index, (chrom, start, end, extra) in enumerate(bed_rows(path)):
            bin_id = extra[0] if extra and extra[0] else "{}:{}-{}".format(chrom, start, end)
            if bin_id in seen:
                raise ValueError("duplicate DELFI bin_id {}".format(bin_id))
            seen.add(bin_id)
            provided_gc = float(extra[1]) if len(extra) > 1 and extra[1] else None
            if provided_gc is not None and not 0 <= provided_gc <= 1:
                raise ValueError("DELFI bin {} GC fraction must be in [0, 1]".format(bin_id))
            row = {"bin_id": bin_id, "chrom": chrom, "start": start, "end": end,
                   "index": index, "provided_gc_fraction": provided_gc}
            self.rows.append(row)
            by_chrom[chrom].append(row)
        self.by_chrom = {}
        self.starts = {}
        for chrom, rows in by_chrom.items():
            rows.sort(key=lambda row: (row["start"], row["end"]))
            for previous, current in zip(rows, rows[1:]):
                if current["start"] < previous["end"]:
                    raise ValueError("overlapping DELFI bins on {} are not supported".format(chrom))
            self.by_chrom[chrom] = rows
            self.starts[chrom] = [row["start"] for row in rows]

    def midpoint_bin(self, chrom, midpoint):
        rows = self.by_chrom.get(chrom, [])
        if not rows:
            return None
        index = bisect.bisect_right(self.starts[chrom], midpoint) - 1
        if index >= 0 and midpoint < rows[index]["end"]:
            return rows[index]
        return None


class GriffinSites:
    def __init__(self, path, window_bp):
        self.rows = []
        by_chrom = defaultdict(list)
        for index, (chrom, start, end, extra) in enumerate(bed_rows(path)):
            if not extra or not extra[0]:
                raise ValueError("Griffin site BED requires TF/site-set name in column 4")
            tf = extra[0]
            site_id = extra[1] if len(extra) > 1 and extra[1] else \
                "{}:{}-{}:{}".format(chrom, start, end, tf)
            row = {"site_id": site_id, "tf": tf, "chrom": chrom,
                   "start": start, "end": end, "center": (start + end) // 2,
                   "index": index}
            self.rows.append(row)
            by_chrom[chrom].append(row)
        self.by_chrom = {}
        self.centers = {}
        self.window_bp = window_bp
        for chrom, rows in by_chrom.items():
            rows.sort(key=lambda row: row["center"])
            self.by_chrom[chrom] = rows
            self.centers[chrom] = [row["center"] for row in rows]

    def candidates(self, chrom, fragment_start, fragment_end):
        rows = self.by_chrom.get(chrom, [])
        if not rows:
            return []
        centers = self.centers[chrom]
        left = bisect.bisect_left(centers, fragment_start - self.window_bp)
        right = bisect.bisect_right(centers, fragment_end + self.window_bp)
        return rows[left:right]


def read_passes(read, options):
    length = int(read.template_length)
    return (read.is_proper_pair and length > 0 and
            read.mapping_quality >= options.min_mapq and
            read.flag & options.exclude_flags == 0 and
            options.min_length <= length <= options.max_length)


def gc_fraction(sequence):
    sequence = sequence.upper()
    valid = sum(sequence.count(base) for base in "ACGT")
    return ((sequence.count("G") + sequence.count("C")) / valid
            if valid else None)


def median(values):
    values = sorted(values)
    if not values:
        return None
    middle = len(values) // 2
    return (values[middle] if len(values) % 2 else
            (values[middle - 1] + values[middle]) / 2)


def mean(values):
    values = [value for value in values if value is not None and finite(value)]
    return sum(values) / len(values) if values else None


def lowess(x, y, span=0.3):
    """Small deterministic local-linear LOWESS without robust iterations."""
    if len(x) != len(y):
        raise ValueError("LOWESS x and y lengths differ")
    n = len(x)
    if n < 3:
        return [None] * n
    neighbours = max(3, min(n, int(math.ceil(span * n))))
    fitted = []
    for target in x:
        distances = sorted(abs(value - target) for value in x)
        bandwidth = distances[neighbours - 1]
        if bandwidth == 0:
            local = [value for value, coordinate in zip(y, x) if coordinate == target]
            fitted.append(mean(local))
            continue
        weights = [(1 - (abs(coordinate - target) / bandwidth) ** 3) ** 3
                   if abs(coordinate - target) <= bandwidth else 0.0
                   for coordinate in x]
        sw = sum(weights)
        sx = sum(weight * coordinate for weight, coordinate in zip(weights, x))
        sy = sum(weight * value for weight, value in zip(weights, y))
        sxx = sum(weight * coordinate * coordinate
                  for weight, coordinate in zip(weights, x))
        sxy = sum(weight * coordinate * value
                  for weight, coordinate, value in zip(weights, x, y))
        denominator = sw * sxx - sx * sx
        if sw == 0:
            fitted.append(None)
        elif abs(denominator) < 1e-12:
            fitted.append(sy / sw)
        else:
            slope = (sw * sxy - sx * sy) / denominator
            intercept = (sy - slope * sx) / sw
            fitted.append(intercept + slope * target)
    return fitted


def load_reference(path, key_fields):
    if not path:
        return {}
    result = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = set(key_fields) | {"reference_mean", "reference_sd"}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise ValueError("reference {} requires columns {}".format(
                path, ", ".join(sorted(required))))
        for row in reader:
            key = tuple(row[field] for field in key_fields)
            if key in result:
                raise ValueError("duplicate reference key {}".format(key))
            reference_mean = float(row["reference_mean"])
            reference_sd = float(row["reference_sd"])
            if not finite(reference_mean) or not finite(reference_sd) or reference_sd <= 0:
                raise ValueError("invalid reference values for {}".format(key))
            result[key] = (reference_mean, reference_sd)
    return result


def delfi_rows(bins, counts, fasta, reference, options):
    rows = []
    eligible = []
    for region in bins.rows:
        short = counts[(region["index"], "short")]
        long = counts[(region["index"], "long")]
        total = short + long
        gc = region["provided_gc_fraction"]
        if gc is None:
            try:
                gc = gc_fraction(fasta.fetch(region["chrom"], region["start"], region["end"]))
            except (KeyError, ValueError, IndexError):
                gc = None
        raw = math.log((short + options.delfi_pseudocount) /
                       (long + options.delfi_pseudocount), 2)
        row = dict(region, sample_id=options.sample_id,
                   assay_partition=options.assay_partition,
                   n_short=short, n_long=long,
                   n_informative=total, gc_fraction=gc, raw_log2_short_long=raw,
                   status="pass" if total >= options.delfi_min_fragments and gc is not None
                   else "insufficient_fragments_or_gc")
        rows.append(row)
        if row["status"] == "pass":
            eligible.append(row)
    if len(eligible) >= options.delfi_min_bins_for_lowess:
        fitted = lowess([row["gc_fraction"] for row in eligible],
                        [row["raw_log2_short_long"] for row in eligible],
                        options.delfi_gc_span)
        centre = median([row["raw_log2_short_long"] for row in eligible])
        for row, expected in zip(eligible, fitted):
            row["gc_expected_log2_ratio"] = expected
            row["gc_corrected_log2_ratio"] = (
                row["raw_log2_short_long"] - expected + centre
                if expected is not None else None)
            ref = reference.get((row["bin_id"],))
            if ref and row["gc_corrected_log2_ratio"] is not None:
                row["reference_mean"], row["reference_sd"] = ref
                row["z_score"] = (row["gc_corrected_log2_ratio"] - ref[0]) / ref[1]
    else:
        for row in eligible:
            row["status"] = "insufficient_bins_for_gc_correction"
    return rows


def griffin_bin_bounds(relative_start, relative_end, window_bp, bin_size):
    n_bins = (2 * window_bp) // bin_size
    left = max(0, int(math.floor((relative_start + window_bp) / bin_size)))
    right = min(n_bins, int(math.ceil((relative_end + window_bp) / bin_size)))
    return left, max(left, right)


def griffin_profiles(sites, differences, fasta, reference, options):
    n_bins = (2 * options.griffin_window_bp) // options.griffin_bin_size_bp
    site_counts = Counter(row["tf"] for row in sites.rows)
    gc_sums = defaultdict(lambda: [0.0] * n_bins)
    gc_counts = defaultdict(lambda: [0] * n_bins)
    for site in sites.rows:
        window_start = site["center"] - options.griffin_window_bp
        window_end = site["center"] + options.griffin_window_bp
        if window_start < 0:
            continue
        try:
            sequence = fasta.fetch(site["chrom"], window_start, window_end)
        except (KeyError, ValueError, IndexError):
            continue
        if len(sequence) != 2 * options.griffin_window_bp:
            continue
        for index in range(n_bins):
            start = index * options.griffin_bin_size_bp
            gc = gc_fraction(sequence[start:start + options.griffin_bin_size_bp])
            if gc is not None:
                gc_sums[site["tf"]][index] += gc
                gc_counts[site["tf"]][index] += 1

    profiles = []
    features = []
    for tf in sorted(site_counts):
        running = 0.0
        coverage = []
        for delta in differences[tf][:-1]:
            running += delta
            coverage.append(running / site_counts[tf])
        gc_values = [(gc_sums[tf][index] / gc_counts[tf][index]
                      if gc_counts[tf][index] else None) for index in range(n_bins)]
        valid = [index for index in range(n_bins) if gc_values[index] is not None]
        expected_by_index = {}
        if len(valid) >= options.griffin_min_bins_for_lowess:
            fitted = lowess([gc_values[index] for index in valid],
                            [coverage[index] for index in valid],
                            options.griffin_gc_span)
            expected_by_index.update(zip(valid, fitted))
        scale = mean(coverage)
        corrected = []
        for index, raw in enumerate(coverage):
            expected = expected_by_index.get(index)
            corrected.append(raw / expected * scale
                             if expected is not None and expected > 0 and scale is not None
                             else None)
        corrected_mean = mean(corrected)
        normalised = [(value / corrected_mean if value is not None and corrected_mean else None)
                      for value in corrected]
        for index in range(n_bins):
            profiles.append({
                "sample_id": options.sample_id,
                "assay_partition": options.assay_partition,
                "tf": tf,
                "position_bp": (-options.griffin_window_bp + index *
                                options.griffin_bin_size_bp +
                                options.griffin_bin_size_bp / 2),
                "n_sites": site_counts[tf],
                "raw_coverage_per_site": coverage[index],
                "mean_site_gc": gc_values[index],
                "gc_expected_coverage": expected_by_index.get(index),
                "gc_corrected_coverage": corrected[index],
                "normalised_coverage": normalised[index],
            })
        usable = [(profiles[-n_bins + index]["position_bp"], normalised[index])
                  for index in range(n_bins) if normalised[index] is not None]
        central = [value for position, value in usable
                   if abs(position) <= options.griffin_central_bp]
        amplitude_region = [value for position, value in usable
                            if abs(position) <= options.griffin_amplitude_bp]
        feature_values = {
            "central_coverage": mean(central),
            "mean_coverage": mean(coverage),
            "amplitude": (max(amplitude_region) - min(amplitude_region)
                          if amplitude_region else None),
        }
        for feature, value in feature_values.items():
            status = ("insufficient_sites" if site_counts[tf] < options.griffin_min_sites else
                      "pass" if value is not None else "gc_correction_failed")
            row = {"sample_id": options.sample_id,
                   "assay_partition": options.assay_partition,
                   "tf": tf, "feature": feature,
                   "value": value, "n_sites": site_counts[tf],
                   "status": status}
            ref = reference.get((tf, feature))
            if ref and value is not None and status == "pass":
                row["reference_mean"], row["reference_sd"] = ref
                row["z_score"] = (value - ref[0]) / ref[1]
            features.append(row)
    return profiles, features


def summarise(delfi, griffin, options):
    z_delfi = [row.get("z_score") for row in delfi if finite(row.get("z_score"))]
    z_griffin = [row.get("z_score") for row in griffin if finite(row.get("z_score"))]
    return {
        "sample_id": options.sample_id,
        "assay_partition": options.assay_partition,
        "delfi_status": ("not_requested" if not options.enable_delfi else
                         "pass" if z_delfi else
                         "reference_not_supplied_or_no_usable_bins"),
        "delfi_n_bins": len(delfi),
        "delfi_n_zscored_bins": len(z_delfi),
        "delfi_mean_abs_z": mean([abs(value) for value in z_delfi]),
        "delfi_rms_z": (math.sqrt(mean([value * value for value in z_delfi]))
                        if z_delfi else None),
        "griffin_status": ("not_requested" if not options.enable_griffin else
                           "pass" if z_griffin else
                           "reference_not_supplied_or_no_usable_features"),
        "griffin_n_features": len(griffin),
        "griffin_n_zscored_features": len(z_griffin),
        "griffin_mean_abs_z": mean([abs(value) for value in z_griffin]),
    }


DELFI_FIELDS = [
    "sample_id", "assay_partition", "bin_id", "chrom", "start", "end", "n_short", "n_long",
    "n_informative", "gc_fraction", "raw_log2_short_long",
    "gc_expected_log2_ratio", "gc_corrected_log2_ratio", "reference_mean",
    "reference_sd", "z_score", "status"
]
GRIFFIN_PROFILE_FIELDS = [
    "sample_id", "assay_partition", "tf", "position_bp", "n_sites", "raw_coverage_per_site",
    "mean_site_gc", "gc_expected_coverage", "gc_corrected_coverage",
    "normalised_coverage"
]
GRIFFIN_FEATURE_FIELDS = [
    "sample_id", "assay_partition", "tf", "feature", "value", "n_sites", "reference_mean",
    "reference_sd", "z_score", "status"
]
REGIONAL_SUMMARY_FIELDS = [
    "sample_id", "assay_partition", "delfi_status", "delfi_n_bins", "delfi_n_zscored_bins",
    "delfi_mean_abs_z", "delfi_rms_z", "griffin_status",
    "griffin_n_features", "griffin_n_zscored_features", "griffin_mean_abs_z"
]


def run(options):
    if not options.enable_delfi and not options.enable_griffin:
        write_tsv(options.delfi_output, [], DELFI_FIELDS)
        write_tsv(options.griffin_profile_output, [], GRIFFIN_PROFILE_FIELDS)
        write_tsv(options.griffin_features_output, [], GRIFFIN_FEATURE_FIELDS)
        write_tsv(options.summary_output, [summarise([], [], options)],
                  REGIONAL_SUMMARY_FIELDS)
        return

    import pysam

    blacklist = Blacklist(options.blacklist)
    bins = DelfiBins(options.delfi_bins) if options.enable_delfi else None
    sites = GriffinSites(options.griffin_sites, options.griffin_window_bp) \
        if options.enable_griffin else None
    delfi_counts = Counter()
    n_griffin_bins = ((2 * options.griffin_window_bp) //
                      options.griffin_bin_size_bp)
    differences = defaultdict(lambda: [0.0] * (n_griffin_bins + 1))
    with pysam.AlignmentFile(options.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if not read_passes(read, options):
                continue
            start = int(read.reference_start)
            end = start + int(read.template_length)
            chrom = read.reference_name
            if blacklist.overlaps(chrom, start, end):
                continue
            length = int(read.template_length)
            if bins:
                region = bins.midpoint_bin(chrom, (start + end) // 2)
                if region:
                    if options.short_min <= length <= options.short_max:
                        delfi_counts[(region["index"], "short")] += 1
                    elif options.long_min <= length <= options.long_max:
                        delfi_counts[(region["index"], "long")] += 1
            if sites:
                for site in sites.candidates(chrom, start, end):
                    relative_start = max(-options.griffin_window_bp, start - site["center"])
                    relative_end = min(options.griffin_window_bp, end - site["center"])
                    left, right = griffin_bin_bounds(
                        relative_start, relative_end, options.griffin_window_bp,
                        options.griffin_bin_size_bp)
                    if right > left:
                        differences[site["tf"]][left] += 1.0
                        differences[site["tf"]][right] -= 1.0

    with pysam.FastaFile(options.fasta) as fasta:
        delfi = (delfi_rows(
            bins, delfi_counts, fasta,
            load_reference(options.delfi_reference, ["bin_id"]), options)
            if bins else [])
        griffin_profiles_rows, griffin_features_rows = (griffin_profiles(
            sites, differences, fasta,
            load_reference(options.griffin_reference, ["tf", "feature"]), options)
            if sites else ([], []))
    summary = summarise(delfi, griffin_features_rows, options)
    write_tsv(options.delfi_output, delfi, DELFI_FIELDS)
    write_tsv(options.griffin_profile_output, griffin_profiles_rows,
              GRIFFIN_PROFILE_FIELDS)
    write_tsv(options.griffin_features_output, griffin_features_rows,
              GRIFFIN_FEATURE_FIELDS)
    write_tsv(options.summary_output, [summary], REGIONAL_SUMMARY_FIELDS)


def parser():
    result = argparse.ArgumentParser()
    result.add_argument("--sample-id", required=True)
    result.add_argument("--assay-partition", required=True)
    result.add_argument("--bam", required=True)
    result.add_argument("--fasta", required=True)
    result.add_argument("--blacklist", required=True)
    result.add_argument("--enable-delfi", type=parse_bool, default=False)
    result.add_argument("--delfi-bins", default="")
    result.add_argument("--delfi-reference", default="")
    result.add_argument("--enable-griffin", type=parse_bool, default=False)
    result.add_argument("--griffin-sites", default="")
    result.add_argument("--griffin-reference", default="")
    result.add_argument("--min-mapq", type=int, default=30)
    result.add_argument("--exclude-flags", type=int, default=3852)
    result.add_argument("--min-length", type=int, default=50)
    result.add_argument("--max-length", type=int, default=250)
    result.add_argument("--short-min", type=int, default=100)
    result.add_argument("--short-max", type=int, default=150)
    result.add_argument("--long-min", type=int, default=151)
    result.add_argument("--long-max", type=int, default=220)
    result.add_argument("--delfi-pseudocount", type=float, default=0.5)
    result.add_argument("--delfi-min-fragments", type=int, default=20)
    result.add_argument("--delfi-min-bins-for-lowess", type=int, default=20)
    result.add_argument("--delfi-gc-span", type=float, default=0.3)
    result.add_argument("--griffin-window-bp", type=int, default=1000)
    result.add_argument("--griffin-bin-size-bp", type=int, default=10)
    result.add_argument("--griffin-central-bp", type=int, default=30)
    result.add_argument("--griffin-amplitude-bp", type=int, default=250)
    result.add_argument("--griffin-min-bins-for-lowess", type=int, default=20)
    result.add_argument("--griffin-min-sites", type=int, default=100)
    result.add_argument("--griffin-gc-span", type=float, default=0.3)
    result.add_argument("--delfi-output", required=True)
    result.add_argument("--griffin-profile-output", required=True)
    result.add_argument("--griffin-features-output", required=True)
    result.add_argument("--summary-output", required=True)
    return result


def main(argv=None):
    options = parser().parse_args(argv)
    if options.enable_delfi and not options.delfi_bins:
        raise ValueError("enable-delfi requires --delfi-bins")
    if options.enable_griffin and not options.griffin_sites:
        raise ValueError("enable-griffin requires --griffin-sites")
    if not 0 < options.delfi_gc_span <= 1 or not 0 < options.griffin_gc_span <= 1:
        raise ValueError("GC LOWESS spans must be in (0, 1]")
    if options.griffin_bin_size_bp <= 0 or options.griffin_window_bp <= 0:
        raise ValueError("Griffin window and bin size must be positive")
    if (2 * options.griffin_window_bp) % options.griffin_bin_size_bp:
        raise ValueError("2 * Griffin window must be divisible by bin size")
    if options.min_length > options.max_length:
        raise ValueError("min-length cannot exceed max-length")
    if options.delfi_pseudocount <= 0:
        raise ValueError("DELFI pseudocount must be positive")
    if options.delfi_min_fragments < 0 or options.delfi_min_bins_for_lowess < 3:
        raise ValueError("DELFI minimum fragments/bins are invalid")
    if options.griffin_min_sites < 1 or options.griffin_min_bins_for_lowess < 3:
        raise ValueError("Griffin minimum sites/bins are invalid")
    if options.griffin_central_bp > options.griffin_window_bp or options.griffin_amplitude_bp > options.griffin_window_bp:
        raise ValueError("Griffin feature windows cannot exceed the profile window")
    run(options)


if __name__ == "__main__":
    main()
