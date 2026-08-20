#!/usr/bin/env python3

"""Build DELFI and Griffin reference tables from independent control outputs."""

import argparse
import csv
import math
import os
from collections import defaultdict


def finite(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def read_values(paths, key_fields, value_field):
    values = defaultdict(dict)
    for path in paths:
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = set(key_fields) | {"sample_id", value_field}
            if not reader.fieldnames or not required <= set(reader.fieldnames):
                raise ValueError("{} requires columns {}".format(
                    path, ", ".join(sorted(required))))
            for row in reader:
                if "status" in row and row["status"] != "pass":
                    continue
                if not finite(row[value_field]):
                    continue
                key = tuple(row[field] for field in key_fields)
                sample = row["sample_id"]
                if sample in values[key]:
                    raise ValueError("duplicate {} value for sample {}".format(key, sample))
                values[key][sample] = float(row[value_field])
    return values


def reference_rows(values, key_fields, min_samples):
    result = []
    for key in sorted(values):
        observations = list(values[key].values())
        if len(observations) < min_samples:
            continue
        centre = sum(observations) / len(observations)
        variance = sum((value - centre) ** 2 for value in observations) / (len(observations) - 1)
        sd = math.sqrt(variance)
        if sd <= 0:
            continue
        row = dict(zip(key_fields, key))
        row.update(reference_mean=centre, reference_sd=sd,
                   n_reference_samples=len(observations))
        result.append(row)
    return result


def write(path, rows, fields):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields,
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def comma_paths(value):
    return [path for path in value.split(",") if path]


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--delfi-inputs", default="")
    parser.add_argument("--griffin-inputs", default="")
    parser.add_argument("--delfi-output", default="")
    parser.add_argument("--griffin-output", default="")
    parser.add_argument("--min-samples", type=int, default=10)
    options = parser.parse_args(argv)
    if options.min_samples < 2:
        raise ValueError("min-samples must be at least two")
    if options.delfi_inputs:
        if not options.delfi_output:
            raise ValueError("DELFI inputs require --delfi-output")
        values = read_values(comma_paths(options.delfi_inputs), ["bin_id"],
                             "gc_corrected_log2_ratio")
        rows = reference_rows(values, ["bin_id"], options.min_samples)
        write(options.delfi_output, rows,
              ["bin_id", "reference_mean", "reference_sd", "n_reference_samples"])
    if options.griffin_inputs:
        if not options.griffin_output:
            raise ValueError("Griffin inputs require --griffin-output")
        values = read_values(comma_paths(options.griffin_inputs), ["tf", "feature"],
                             "value")
        rows = reference_rows(values, ["tf", "feature"], options.min_samples)
        write(options.griffin_output, rows,
              ["tf", "feature", "reference_mean", "reference_sd",
               "n_reference_samples"])


if __name__ == "__main__":
    main()
