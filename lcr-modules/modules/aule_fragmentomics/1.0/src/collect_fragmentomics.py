#!/usr/bin/env python3

import argparse
import bisect
import csv
import gzip
import itertools
import math
import os
from collections import Counter, defaultdict


DNA = "ACGT"
MOTIFS = ["".join(x) for x in itertools.product(DNA, repeat=4)]


def parse_bool(value):
    value = str(value).strip().lower()
    if value in {"true", "1", "yes"}:
        return True
    if value in {"false", "0", "no"}:
        return False
    raise argparse.ArgumentTypeError("expected true or false")


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def finite(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def display(value):
    if value is None or not finite(value) and not isinstance(value, (str, bool)):
        return "NA"
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    return value


def ensure_parent(path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)


def write_tsv(path, rows, fieldnames):
    ensure_parent(path)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames,
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: display(row.get(key)) for key in fieldnames})


class IntervalIndex:
    def __init__(self):
        self.intervals = defaultdict(list)
        self.starts = {}

    def add(self, chrom, start, end):
        if end > start:
            self.intervals[chrom].append((start, end))

    def finalise(self):
        for chrom, intervals in list(self.intervals.items()):
            merged = []
            for start, end in sorted(intervals):
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    merged.append((start, end))
            self.intervals[chrom] = merged
            self.starts[chrom] = [start for start, _ in merged]
        return self

    def overlaps(self, chrom, start, end):
        intervals = self.intervals.get(chrom, [])
        if not intervals:
            return False
        index = bisect.bisect_left(self.starts[chrom], end) - 1
        return index >= 0 and intervals[index][1] > start


def read_blacklist(path):
    index = IntervalIndex()
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError("blacklist line {} has fewer than three columns".format(line_number))
            index.add(fields[0], int(fields[1]), int(fields[2]))
    return index.finalise()


def read_passes(read, min_mapq, exclude_flags, min_length, max_length):
    length = int(read.template_length)
    return (read.is_proper_pair and length > 0 and
            read.mapping_quality >= min_mapq and
            read.flag & exclude_flags == 0 and
            min_length <= length <= max_length)


def allele_read_passes(read, min_mapq, exclude_flags, min_length, max_length):
    length = abs(int(read.template_length))
    return (read.is_proper_pair and length > 0 and
            read.mapping_quality >= min_mapq and
            read.flag & exclude_flags == 0 and
            min_length <= length <= max_length)


def fragment_bounds(read):
    length = int(read.template_length)
    if length > 0:
        start = int(read.reference_start)
        return start, start + length
    end = int(read.reference_end)
    return end + length, end


def fragment_end_motifs(fasta, chrom, start, end):
    if start < 0 or end - start < 4:
        return []
    try:
        left = fasta.fetch(chrom, start, start + 4).upper()
        right = reverse_complement(fasta.fetch(chrom, end - 4, end).upper())
    except (KeyError, ValueError, IndexError):
        return []
    return [motif for motif in (left, right)
            if len(motif) == 4 and set(motif) <= set(DNA)]


def weighted_mean(histogram):
    total = sum(histogram.values())
    return sum(length * count for length, count in histogram.items()) / total if total else None


def weighted_median(histogram):
    total = sum(histogram.values())
    if not total:
        return None
    left_rank = (total - 1) // 2
    right_rank = total // 2
    cumulative = 0
    values = []
    for length, count in sorted(histogram.items()):
        previous = cumulative
        cumulative += count
        if previous <= left_rank < cumulative:
            values.append(length)
        if previous <= right_rank < cumulative:
            values.append(length)
        if len(values) == 2:
            return sum(values) / 2
    return None


def pearson(x, y):
    if len(x) != len(y) or len(x) < 3:
        return None
    mean_x = sum(x) / len(x)
    mean_y = sum(y) / len(y)
    dx = [value - mean_x for value in x]
    dy = [value - mean_y for value in y]
    denom = math.sqrt(sum(value * value for value in dx) *
                      sum(value * value for value in dy))
    return sum(a * b for a, b in zip(dx, dy)) / denom if denom else None


def periodicity_score(histogram, min_length, max_length, min_period, max_period):
    values = [float(histogram.get(length, 0))
              for length in range(min_length, max_length + 1)]
    if len(values) < 2 * max_period + 3 or sum(values) == 0:
        return None, None
    radius = 10
    smooth = []
    for index in range(len(values)):
        lo = max(0, index - radius)
        hi = min(len(values), index + radius + 1)
        smooth.append(sum(values[lo:hi]) / (hi - lo))
    residual = [value - trend for value, trend in zip(values, smooth)]
    candidates = []
    for lag in range(min_period, max_period + 1):
        score = pearson(residual[:-lag], residual[lag:])
        if score is not None:
            candidates.append((score, lag))
    return max(candidates) if candidates else (None, None)


def normalised_entropy(counts, categories):
    total = sum(counts.values())
    if total <= 0 or categories <= 1:
        return None
    entropy = -sum((count / total) * math.log(count / total)
                   for count in counts.values() if count > 0)
    return entropy / math.log(categories)


def size_metrics(histogram, options):
    total = sum(histogram.values())
    short = sum(count for length, count in histogram.items()
                if options.short_min <= length <= options.short_max)
    long = sum(count for length, count in histogram.items()
               if options.long_min <= length <= options.long_max)
    mode = min((length for length, count in histogram.items()
                if count == max(histogram.values())), default=None) if total else None
    amplitude, period = periodicity_score(
        histogram, options.min_length,
        min(options.max_length, options.periodicity_max_length),
        options.periodicity_min_bp, options.periodicity_max_bp)
    return {
        "n_fragments": total,
        "mean_fragment_length": weighted_mean(histogram),
        "median_fragment_length": weighted_median(histogram),
        "modal_fragment_length": mode,
        "n_short_100_150": short,
        "n_long_151_220": long,
        "p_short": short / total if total else None,
        "short_long_ratio": short / long if long else None,
        "periodicity_amplitude": amplitude,
        "periodicity_bp": period,
    }


def collect_global(bam_path, fasta_path, blacklist, options):
    import pysam

    lengths = Counter()
    motifs = Counter()
    diagnostics = Counter()
    with pysam.AlignmentFile(bam_path, "rb") as bam, pysam.FastaFile(fasta_path) as fasta:
        for read in bam.fetch(until_eof=True):
            diagnostics["alignments_seen"] += 1
            if not read_passes(read, options.min_mapq, options.exclude_flags,
                               options.min_length, options.max_length):
                diagnostics["alignments_filtered"] += 1
                continue
            start, end = fragment_bounds(read)
            chrom = read.reference_name
            if blacklist.overlaps(chrom, start, end):
                diagnostics["blacklisted_fragments"] += 1
                continue
            lengths[int(read.template_length)] += 1
            end_motifs = fragment_end_motifs(fasta, chrom, start, end)
            motifs.update(end_motifs)
            diagnostics["motif_ends_missing"] += 2 - len(end_motifs)
    return lengths, motifs, diagnostics


def base_at_reference_position(read, position, min_baseq=0):
    sequence = read.query_sequence
    if not sequence:
        return None
    for query_position, reference_position in read.get_aligned_pairs(matches_only=False):
        if reference_position == position:
            if query_position is None:
                return None
            qualities = read.query_qualities
            if (min_baseq > 0 and
                    (qualities is None or query_position >= len(qualities) or
                     qualities[query_position] < min_baseq)):
                return None
            return sequence[query_position].upper()
    return None


def ks_2sample(x, y):
    if not x or not y:
        return None, None
    x = sorted(x)
    y = sorted(y)
    values = sorted(set(x + y))
    distance = max(abs(bisect.bisect_right(x, value) / len(x) -
                       bisect.bisect_right(y, value) / len(y)) for value in values)
    if distance == 0:
        return 0.0, 1.0
    effective = len(x) * len(y) / (len(x) + len(y))
    corrected = (math.sqrt(effective) + 0.12 + 0.11 / math.sqrt(effective)) * distance
    p_value = 2 * sum(((-1) ** (term - 1)) * math.exp(-2 * term * term * corrected * corrected)
                      for term in range(1, 101))
    return distance, min(1.0, max(0.0, p_value))


def distribution_metrics(mutant, wild_type):
    mutant_hist = Counter(mutant)
    wild_hist = Counter(wild_type)
    ks_distance, ks_p = ks_2sample(mutant, wild_type)
    mutant_mean = weighted_mean(mutant_hist)
    wild_mean = weighted_mean(wild_hist)
    return {
        "n_mutant": len(mutant),
        "n_wild_type": len(wild_type),
        "mutant_mean_length": mutant_mean,
        "wild_type_mean_length": wild_mean,
        "mutant_median_length": weighted_median(mutant_hist),
        "wild_type_median_length": weighted_median(wild_hist),
        "mean_length_shift_mutant_minus_wt": (
            mutant_mean - wild_mean if mutant_mean is not None and wild_mean is not None else None),
        "ks_distance": ks_distance,
        "ks_p_value": ks_p,
    }


VARIANT_FIELDS = [
    "sample_id", "variant_id", "chrom", "pos_1based", "ref", "alt", "status",
    "n_mutant", "n_wild_type", "mutant_mean_length", "wild_type_mean_length",
    "mutant_median_length", "wild_type_median_length",
    "mean_length_shift_mutant_minus_wt", "ks_distance", "ks_p_value"
]


def collect_alleles(bam_path, vcf_path, blacklist, options):
    import pysam

    rows = []
    fragment_calls = defaultdict(list)
    conflicting_fragment_names = set()
    supported_variants = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam, pysam.VariantFile(vcf_path) as vcf:
        for record in vcf:
            variant_id = record.id or "{}:{}:{}>{}".format(
                record.chrom, record.pos, record.ref, ",".join(record.alts or []))
            base_row = {
                "sample_id": options.sample_id,
                "variant_id": variant_id,
                "chrom": record.chrom,
                "pos_1based": record.pos,
                "ref": record.ref,
                "alt": ",".join(record.alts or []),
            }
            filters = set(record.filter.keys())
            if options.pass_variants_only and filters and "PASS" not in filters:
                rows.append(dict(base_row, status="filtered_variant"))
                continue
            if not record.alts or len(record.alts) != 1 or len(record.ref) != 1 or len(record.alts[0]) != 1:
                rows.append(dict(base_row, status="unsupported_non_biallelic_snv"))
                continue
            if blacklist.overlaps(record.chrom, record.start, record.stop):
                rows.append(dict(base_row, status="blacklisted_variant"))
                continue
            mutant = []
            wild_type = []
            variant_fragment_calls = defaultdict(list)
            try:
                reads = bam.fetch(record.chrom, record.start, record.stop)
                for read in reads:
                    if not allele_read_passes(
                            read, options.min_mapq, options.exclude_flags,
                            options.min_length, options.max_length):
                        continue
                    start, end = fragment_bounds(read)
                    if blacklist.overlaps(record.chrom, start, end):
                        continue
                    base = base_at_reference_position(
                        read, record.start, options.min_baseq)
                    if base == record.alts[0].upper():
                        allele = "mutant"
                    elif base == record.ref.upper():
                        allele = "wild_type"
                    else:
                        continue
                    variant_fragment_calls[read.query_name].append(
                        (allele, abs(int(read.template_length))))
            except ValueError:
                rows.append(dict(base_row, status="contig_missing_from_bam"))
                continue
            for query_name, calls in variant_fragment_calls.items():
                alleles = {allele for allele, _ in calls}
                lengths = {length for _, length in calls}
                if len(alleles) != 1 or len(lengths) != 1:
                    conflicting_fragment_names.add(query_name)
                    continue
                allele = next(iter(alleles))
                length = next(iter(lengths))
                (mutant if allele == "mutant" else wild_type).append(length)
                fragment_calls[query_name].append((allele, length))
            supported_variants += 1
            metrics = distribution_metrics(mutant, wild_type)
            status = "pass" if mutant and wild_type else "insufficient_allele_observations"
            rows.append(dict(base_row, status=status, **metrics))

    mutant_unique = []
    wild_unique = []
    for query_name, calls in fragment_calls.items():
        alleles = {allele for allele, _ in calls}
        lengths = {length for _, length in calls}
        if len(alleles) != 1 or len(lengths) != 1:
            conflicting_fragment_names.add(query_name)
            continue
        allele = next(iter(alleles))
        length = next(iter(lengths))
        (mutant_unique if allele == "mutant" else wild_unique).append(length)
    summary = distribution_metrics(mutant_unique, wild_unique)
    summary.update({
        "n_supported_variants": supported_variants,
        "n_conflicting_fragments": len(conflicting_fragment_names),
        "raw_variant_vaf": (len(mutant_unique) / (len(mutant_unique) + len(wild_unique))
                            if mutant_unique or wild_unique else None),
    })
    if supported_variants == 0:
        summary["allele_status"] = "no_supported_variants"
    elif len(mutant_unique) < options.min_mutant_fragments:
        summary["allele_status"] = "insufficient_mutant_fragments"
    elif not wild_unique:
        summary["allele_status"] = "no_wild_type_fragments"
    else:
        summary["allele_status"] = "pass"
    return summary, rows, Counter(mutant_unique), Counter(wild_unique)


def load_size_weights(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "length" not in reader.fieldnames:
            raise ValueError("size-weight table requires a length column")
        weight_column = "likelihood_ratio" if "likelihood_ratio" in reader.fieldnames else "weight"
        if weight_column not in reader.fieldnames:
            raise ValueError("size-weight table requires likelihood_ratio or weight")
        weights = {}
        for row in reader:
            length = int(row["length"])
            weight = float(row[weight_column])
            if not math.isfinite(weight) or weight < 0:
                raise ValueError("invalid weight for fragment length {}".format(length))
            if length in weights and weights[length] != weight:
                raise ValueError("duplicate conflicting weight for length {}".format(length))
            weights[length] = weight
    return weights


def apply_weights(mutant_hist, wild_hist, weights):
    missing = sorted({length for length in mutant_hist.keys() | wild_hist.keys()
                      if length not in weights})
    if missing:
        return {
            "weighting_status": "incomplete_weight_table",
            "n_weight_lengths_missing": len(missing),
            "weighted_mutant_count": None,
            "weighted_wild_type_count": None,
            "weighted_variant_vaf": None,
        }
    mutant = sum(count * weights[length] for length, count in mutant_hist.items())
    wild = sum(count * weights[length] for length, count in wild_hist.items())
    return {
        "weighting_status": "pass",
        "n_weight_lengths_missing": 0,
        "weighted_mutant_count": mutant,
        "weighted_wild_type_count": wild,
        "weighted_variant_vaf": mutant / (mutant + wild) if mutant + wild else None,
    }


SUMMARY_FIELDS = [
    "sample_id", "assay_partition", "n_alignments_seen", "n_alignments_filtered",
    "n_blacklisted_fragments", "n_fragments", "mean_fragment_length",
    "median_fragment_length", "modal_fragment_length", "n_short_100_150",
    "n_long_151_220", "p_short", "short_long_ratio", "periodicity_amplitude",
    "periodicity_bp", "n_motif_ends", "n_motif_ends_missing",
    "motif_diversity", "allele_status", "n_supported_variants",
    "n_conflicting_fragments", "n_mutant", "n_wild_type",
    "mutant_mean_length", "wild_type_mean_length", "mutant_median_length",
    "wild_type_median_length", "mean_length_shift_mutant_minus_wt",
    "ks_distance", "ks_p_value", "raw_variant_vaf", "weighting_status",
    "n_weight_lengths_missing", "weighted_mutant_count",
    "weighted_wild_type_count", "weighted_variant_vaf", "min_mapq", "min_baseq",
    "exclude_flags", "min_fragment_length", "max_fragment_length",
    "blacklist_bed"
]


def run(options):
    blacklist = read_blacklist(options.blacklist)
    size_hist, motif_counts, diagnostics = collect_global(
        options.bam, options.fasta, blacklist, options)
    summary = {"sample_id": options.sample_id,
               "assay_partition": options.assay_partition}
    summary.update(size_metrics(size_hist, options))
    summary.update({
        "n_alignments_seen": diagnostics["alignments_seen"],
        "n_alignments_filtered": diagnostics["alignments_filtered"],
        "n_blacklisted_fragments": diagnostics["blacklisted_fragments"],
        "n_motif_ends": sum(motif_counts.values()),
        "n_motif_ends_missing": diagnostics["motif_ends_missing"],
        "motif_diversity": normalised_entropy(motif_counts, 256),
        "min_mapq": options.min_mapq,
        "min_baseq": options.min_baseq,
        "exclude_flags": options.exclude_flags,
        "min_fragment_length": options.min_length,
        "max_fragment_length": options.max_length,
        "blacklist_bed": options.blacklist,
    })

    if options.variant_vcf:
        allele, variant_rows, mutant_hist, wild_hist = collect_alleles(
            options.bam, options.variant_vcf, blacklist, options)
        summary.update(allele)
    else:
        variant_rows = []
        mutant_hist = Counter()
        wild_hist = Counter()
        summary["allele_status"] = "not_requested"

    if options.size_weights:
        summary.update(apply_weights(
            mutant_hist, wild_hist, load_size_weights(options.size_weights)))
    else:
        summary["weighting_status"] = "not_requested"

    write_tsv(options.summary, [summary], SUMMARY_FIELDS)
    total = sum(size_hist.values())
    write_tsv(options.size_hist, [
        {"sample_id": options.sample_id, "length": length,
         "count": size_hist.get(length, 0),
         "frequency": size_hist.get(length, 0) / total if total else None}
        for length in range(options.min_length, options.max_length + 1)
    ], ["sample_id", "length", "count", "frequency"])
    motif_total = sum(motif_counts.values())
    write_tsv(options.motifs, [
        {"sample_id": options.sample_id, "motif": motif,
         "count": motif_counts.get(motif, 0),
         "frequency": motif_counts.get(motif, 0) / motif_total if motif_total else None}
        for motif in MOTIFS
    ], ["sample_id", "motif", "count", "frequency"])
    allele_rows = []
    for allele, histogram in (("mutant", mutant_hist), ("wild_type", wild_hist)):
        allele_total = sum(histogram.values())
        for length in range(options.min_length, options.max_length + 1):
            allele_rows.append({
                "sample_id": options.sample_id, "allele": allele, "length": length,
                "count": histogram.get(length, 0),
                "frequency": histogram.get(length, 0) / allele_total if allele_total else None,
            })
    write_tsv(options.allele_hist, allele_rows,
              ["sample_id", "allele", "length", "count", "frequency"])
    write_tsv(options.variants, variant_rows, VARIANT_FIELDS)


def parser():
    result = argparse.ArgumentParser()
    result.add_argument("--sample-id", required=True)
    result.add_argument("--assay-partition", required=True)
    result.add_argument("--bam", required=True)
    result.add_argument("--fasta", required=True)
    result.add_argument("--blacklist", required=True)
    result.add_argument("--variant-vcf", default="")
    result.add_argument("--size-weights", default="")
    result.add_argument("--pass-variants-only", type=parse_bool, default=True)
    result.add_argument("--min-mapq", type=int, default=30)
    result.add_argument("--min-baseq", type=int, default=30)
    result.add_argument("--exclude-flags", type=int, default=3852)
    result.add_argument("--min-length", type=int, default=50)
    result.add_argument("--max-length", type=int, default=250)
    result.add_argument("--short-min", type=int, default=100)
    result.add_argument("--short-max", type=int, default=150)
    result.add_argument("--long-min", type=int, default=151)
    result.add_argument("--long-max", type=int, default=220)
    result.add_argument("--periodicity-max-length", type=int, default=150)
    result.add_argument("--periodicity-min-bp", type=int, default=8)
    result.add_argument("--periodicity-max-bp", type=int, default=12)
    result.add_argument("--min-mutant-fragments", type=int, default=5)
    result.add_argument("--summary", required=True)
    result.add_argument("--size-hist", required=True)
    result.add_argument("--motifs", required=True)
    result.add_argument("--allele-hist", required=True)
    result.add_argument("--variants", required=True)
    return result


def main(argv=None):
    options = parser().parse_args(argv)
    if options.min_length > options.max_length:
        raise ValueError("min-length cannot exceed max-length")
    if options.periodicity_min_bp > options.periodicity_max_bp:
        raise ValueError("periodicity-min-bp cannot exceed periodicity-max-bp")
    if options.size_weights and not options.variant_vcf:
        raise ValueError("size weights require a variant VCF")
    run(options)


if __name__ == "__main__":
    main()
