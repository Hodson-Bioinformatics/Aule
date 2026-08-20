#!/usr/bin/env python3

import argparse
import csv
import os


def read_regional(paths):
    result = {}
    fields = []
    for path in paths:
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or "sample_id" not in reader.fieldnames:
                raise ValueError("regional summary {} requires sample_id".format(path))
            for field in reader.fieldnames:
                if field != "sample_id" and field not in fields:
                    fields.append(field)
            rows = list(reader)
            if len(rows) != 1:
                raise ValueError("regional summary {} must contain one row".format(path))
            sample_id = rows[0]["sample_id"]
            if sample_id in result:
                raise ValueError("duplicate regional sample_id {}".format(sample_id))
            result[sample_id] = rows[0]
    return result, fields


def aggregate(paths, regional_paths=None):
    rows = []
    fields = []
    for path in paths:
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                raise ValueError("missing header in {}".format(path))
            for field in reader.fieldnames:
                if field not in fields:
                    fields.append(field)
            rows.extend(reader)
    if "sample_id" not in fields:
        raise ValueError("per-sample summaries require sample_id")
    sample_ids = [row.get("sample_id", "") for row in rows]
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("duplicate sample_id in cohort summaries")
    if regional_paths:
        regional, regional_fields = read_regional(regional_paths)
        missing = sorted(set(sample_ids) - set(regional))
        unexpected = sorted(set(regional) - set(sample_ids))
        if missing or unexpected:
            raise ValueError("regional/global sample mismatch: missing={} unexpected={}".format(
                ",".join(missing), ",".join(unexpected)))
        fields.extend(field for field in regional_fields if field not in fields)
        for row in rows:
            row.update({field: regional[row["sample_id"]].get(field, "")
                        for field in regional_fields})
    return sorted(rows, key=lambda row: row["sample_id"]), fields


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", required=True)
    parser.add_argument("--regional-inputs", default="")
    parser.add_argument("--output", required=True)
    options = parser.parse_args(argv)
    paths = [path for path in options.inputs.split(",") if path]
    if not paths:
        raise ValueError("no per-sample summaries supplied")
    regional_paths = [path for path in options.regional_inputs.split(",") if path]
    rows, fields = aggregate(paths, regional_paths)
    os.makedirs(os.path.dirname(os.path.abspath(options.output)), exist_ok=True)
    with open(options.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields,
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
