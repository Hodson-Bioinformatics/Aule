#!/usr/bin/env python3

"""
This script creates IGV batch scripts, which automate the
process of taking IGV screenshots for regions of interest
(e.g., point mutations and structural variations). For
more information, check out this link:

    https://software.broadinstitute.org/software/igv/batch.

The following command will print usage information, including
the various command-line options and the list of supported
file formats:

    generate_igv_batch.py --help;

The script takes the following as inputs:

(1) File specifying mutation sites (the list of supported file
    formats is listed in the usage information when you run
    the script with `--help`);
(2) List of BAM files, which can be provided at the command line
    with `--bam_list` or in a table format using `--bam_table`;
(3) (Optional) List of genes to filter on (not all formats are
    supported);
(4) Additional parameters (check out the help page using
    `--help` to find out more).
"""

# Import standard libraries
import os
import warnings
import argparse

# Import third-party libraries
import vcf
import numpy as np
import pandas as pd


# Define list of supported file formats
FILE_FORMATS = ["maf", "bed", "bedpe", "vcf"]


def main():
    """Runs the script."""

    args = parse_arguments()

    args = validate_arguments(args)

    regions = get_regions(
        args.input_file,
        args.file_format,
        padding=args.padding,
        pass_only=args.pass_only,
        gene_list=args.gene_list,
    )

    generate_igv_batch(
        args.file_format,
        regions,
        args.bam_list,
        args.bam_table,
        args.output,
        args.max_height,
        args.snapshot_dir,
        args.genome_build,
    )

    close_files(args)


def parse_arguments():
    """Parses command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "input_file",
        type=argparse.FileType("r"),
        default="-",
        help=f"Input file {FILE_FORMATS}. Can be '-' for stdin.",
    )

    parser.add_argument(
        "--output",
        "-o",
        metavar="OUTPUT_FILE",
        type=argparse.FileType("w"),
        default="-",
        help="Output IGV batch script. Can be '-' for stdout.",
    )

    parser.add_argument(
        "--format",
        "-f",
        choices=FILE_FORMATS,
        dest="file_format",
        help=(
            "Format of the input file. Otherwise, the format "
            "will be inferred from the file extension."
        ),
    )

    parser.add_argument(
        "--bam_table",
        "-bt",
        metavar="BAM_TABLE_FILE",
        type=argparse.FileType("r"),
        help="Two-column tab-delimited table mapping sample IDs to BAM files paths.",
    )

    parser.add_argument(
        "--bam_list",
        "-bl",
        nargs="+",
        metavar=("BAM_FILE_1", "BAM_FILE_2"),
        default=[],
        help="Space-delimited list of BAM files.",
    )

    parser.add_argument(
        "--gene_column",
        "-gt",
        metavar="GENE_COLUMN_FILE",
        type=argparse.FileType("r"),
        help=(
            "Single-column file with gene symbols (no header). "
            "Used to filter input file. Supported formats: MAF."
        ),
    )

    parser.add_argument(
        "--gene_list",
        "-gl",
        nargs="+",
        metavar=("GENE_1", "GENE_2"),
        default=[],
        help=(
            "Space-delimited list of gene symbols. Used to "
            "filter input file. Supported formats: MAF."
        ),
    )

    default_padding = 300
    parser.add_argument(
        "--padding",
        "-p",
        type=int,
        default=default_padding,
        help=(
            "Amount of padding added before and after each locus. "
            f"Default padding is {default_padding}."
        ),
    )

    default_max_height = 400
    parser.add_argument(
        "--max_height",
        "-m",
        type=int,
        default=default_max_height,
        help=(
            "Maximum panel height in IGV. Default max "
            f"height is {default_max_height}."
        ),
    )

    parser.add_argument(
        "--snapshot_dir",
        "-d",
        required=True,
        help=(
            "Directory where all of the IGV snapshots will "
            "be stored once IGV is run subsequently."
        ),
    )

    parser.add_argument(
        "--genome_build",
        "-g",
        required=True,
        help="Specify IGV genome build (e.g., 'hg19' or 'hg38')",
    )

    parser.add_argument(
        "--pass_only",
        "-po",
        action="store_true",
        help=(
            "Filter for records that have passed any filters. "
            "For example, for VCF files, this would filter for "
            "mutations where the FILTER column contains 'PASS'. "
            "Supported formats: VCF."
        ),
    )

    args = parser.parse_args()

    return args


def validate_arguments(args):
    """Validates command-line arguments."""

    os.makedirs(args.snapshot_dir, exist_ok=True)

    if not args.file_format:
        args.file_format = infer_format(args.input_file)

    assert (
        args.file_format in FILE_FORMATS
    ), f"{args.file_format} format is not supported."

    if args.gene_column:
        gene_column_list = read_lines(args.gene_column)
        args.gene_list += gene_column_list

    return args


def infer_format(input_file):
    """Infers file format based on file extension."""

    input_file_name = input_file.name
    _, input_file_ext = os.path.splitext(input_file_name)
    input_file_ext = input_file_ext.lstrip(".")

    if input_file_ext == "gz":
        raise IOError(
            "Cannot handle compressed input files. You may pipe a decompressed "
            "copy to this script while setting the `input` as '-' to read from "
            "stdin. If so, make sure to provide a value for --format."
        )

    assert input_file_ext in FILE_FORMATS, (
        f"Cannot infer file format from file extension ({input_file_ext}). File "
        f"extension must be among {FILE_FORMATS}. Consider providing --format."
    )

    return input_file_ext


def read_lines(input_file):
    """Generates a list of strings from a single-column file (no header)."""
    df = pd.read_table(input_file, header=None, names=["strings"])
    str_list = df.strings.tolist()
    return str_list


def get_regions_vcf(vcf_file, pass_only, padding, **kwargs):
    """Generates a data frame of regions from a VCF file.

    For paired events (like breakends), only one of each pair is returned.
    """

    vcf_reader = vcf.Reader(vcf_file)
    regions = []
    events_seen = set()

    for record in vcf_reader:
        record_regions = []
        template = "{}:{}-{}"

        if pass_only and len(record.FILTER) > 0:
            continue

        if record.ID in events_seen:
            continue

        region1 = template.format(
            record.CHROM, record.POS - padding, record.POS + padding
        )
        record_regions.append(region1)

        if record.is_sv and "END" in record.INFO:
            end = record.INFO["END"][0]
            region2 = template.format(record.CHROM, end - padding, end + padding)
            record_regions.append(region2)

        elif record.is_sv and record.INFO["SVTYPE"] == "BND":
            region2 = template.format(
                record.ALT[0].chr,
                record.ALT[0].pos - padding,
                record.ALT[0].pos + padding,
            )
            record_regions.append(region2)

            # To skip mate event in VCF file
            events_seen.add(record.INFO["MATEID"])

        regions.append(" ".join(record_regions))

    regions_df = pd.DataFrame(regions, columns=["regions"])

    return regions_df


def get_regions_bed(bed_file, padding, **kwargs):
    """Generates a data frame of regions from a BED file."""

    regions = []
    names = []

    for line in bed_file:
        fields = line.rstrip("\n").split("\t")
        chrom, start, end, name = fields[:4]
        start = int(start) - padding
        end = int(end) + padding
        output = f"{chrom}:{start}-{end}"
        regions.append(output)
        names.append(name)

    regions_df = pd.DataFrame({"regions": regions, "region_name": names})

    return regions_df


def get_regions_bedpe(bedpe_file, padding, **kwargs):
    """Generates a data frame of regions from a BEDPE file."""

    regions = []

    for line in bedpe_file:
        if line.startswith("#") or line.startswith("CHROM"):
            continue
        fields = line.rstrip("\n").split("\t")
        chrom1, start1, end1, chrom2, start2, end2 = fields[:6]
        start1 = int(start1) - padding
        end1 = int(end1) + padding
        start2 = int(start2) - padding
        end2 = int(end2) + padding
        region = f"{chrom1}:{start1}-{end1} {chrom2}:{start2}-{end2}"
        regions.append(region)

    regions_df = pd.DataFrame(regions, columns=["regions"])

    return regions_df


def get_regions_maf(maf_file, padding, gene_list, **kwargs):
    """Generates a data frame of regions from a MAF file."""

    maf = pd.read_table(maf_file, comment="#")

    # Check for issues with input file
    columns = [
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Hugo_Symbol",
        "Tumor_Sample_Barcode",
    ]
    for column in columns:
        column_values = maf[column]
        is_any_na = pd.isna(column_values).any()
        assert not is_any_na, (
            f"The '{column}' column contains NA values. This might be caused "
            "by an incorrectly formatted input MAF file. Please ensure that "
            f"all of the following columns have values: {', '.join(columns)}."
            f"Here's a preview of the MAF file after being parsed:\n\n {maf}"
        )

    if gene_list:
        maf = maf[maf.Hugo_Symbol.isin(gene_list)]

    chrom = maf.Chromosome.astype(str)
    start = (maf.Start_Position - padding).astype(str)
    end = (maf.End_Position + padding).astype(str)
    regions = chrom + ":" + start + "-" + end

    regions_df = pd.DataFrame(
        {
            "regions": regions,
            "region_name": maf.Hugo_Symbol,
            "sample_id": maf.Tumor_Sample_Barcode,
        }
    )

    return regions_df


def get_regions(input_file, file_format, **kwargs):
    """Generates a data frame of regions for a number of file formats.

    Returns pandas DataFrame with the following columns:
        1) regions: With the format "{chrom}:{start}-{end}"
        2) region_name: If a name is associated with the region
        3) sample_id: If a sample is associated with the region
    """

    # Associate file formats with function to get regions
    get_regions_functions = {
        "vcf": get_regions_vcf,
        "bed": get_regions_bed,
        "bedpe": get_regions_bedpe,
        "maf": get_regions_maf,
    }

    # Confirm that all supported file formats have a corresponding function
    assert all(
        f in get_regions_functions for f in FILE_FORMATS
    ), "Not all formats in `FILE_FORMATS` implemented yet."

    get_regions = get_regions_functions[file_format]

    return get_regions(input_file, **kwargs)


def output_lines(lines, output):
    """Outputs list of lines to a file in disk."""
    lines.append("")
    text = "\n".join(lines)
    output.write(text)


def generate_igv_batch_header(bam_files, max_height, snapshot_dir, genome_build):
    """Generates a list of lines that form the header of an IGV batch file."""
    lines = []
    for bam_file in bam_files:
        bam_file = os.path.realpath(bam_file)
        lines.append(f"load {bam_file}")
    lines.append(f"maxPanelHeight {max_height}")
    lines.append(f"snapshotDirectory {snapshot_dir}")
    lines.append(f"genome {genome_build}")
    return lines


def generate_igv_batch_per_row(regions, snapshot_filename):
    """Generates a list of lines that form the body of an IGV batch file."""
    lines = []
    lines.append(f"goto {regions}")
    lines.append("sort")
    lines.append("collapse")
    lines.append(f"snapshot {snapshot_filename}")
    lines.append("new")
    return lines


def generate_igv_batch_footer():
    """Generates a list of lines that form the footer of an IGV batch file."""
    lines = []
    lines.append("exit")
    return lines


def generate_bam_dict(bam_table):
    """Maps samples to BAM files from a BAM table."""
    bam_files = pd.read_table(bam_table, names=["sample_id", "bam_file"])
    bam_files = bam_files.drop_duplicates()
    bam_dict = bam_files.groupby("sample_id")["bam_file"].apply(list).to_dict()
    return bam_dict


def generate_igv_batch_per_sample(
    regions, bam_files, output, max_height, snapshot_dir, genome_build
):
    """Generate IGV batch file per sample."""
    lines = []

    if not bam_files:
        return

    for _, row in regions.iterrows():
        filename = []

        filename.append(row.regions)

        header = generate_igv_batch_header(
            bam_files, max_height, snapshot_dir, genome_build
        )
        lines.extend(header)

        if "region_name" in row:
            filename.append(row.region_name)

        filename = "--".join(filename) + ".png"
        filename = filename.replace(" ", "_")

        row_lines = generate_igv_batch_per_row(row.regions, filename)

        lines.extend(row_lines)

    return lines


def generate_igv_batch(
    file_format,
    regions,
    bam_list,
    bam_table,
    output,
    max_height,
    snapshot_dir,
    genome_build,
):
    """Generate IGV batch file."""

    all_lines = []

    if file_format == "maf":

        samples = regions.sample_id.unique()

        assert bam_list or bam_table, (
            "Please provide a list of BAM files using --bam_list (ideal for single-"
            "sample MAF files) or --bam_table (ideal for multi-sample MAF files)."
        )

        if bam_table:
            bam_dict = generate_bam_dict(bam_table)
        elif len(samples) > 1:
            warnings.warn(
                "Consider using --bam_table for MAF files with multiple samples."
            )

        for sample in samples:

            sample_regions = regions[regions.sample_id == sample]

            bam_files = bam_list
            if bam_table:
                bam_files.extend(bam_dict.get(sample, []))

            sample_snapshot_dir = os.path.join(snapshot_dir, sample, "")

            lines = generate_igv_batch_per_sample(
                sample_regions,
                bam_files,
                output,
                max_height,
                sample_snapshot_dir,
                genome_build,
            )

            all_lines.extend(lines)

    else:

        assert bam_list, "Please provide list of BAM files using --bam_list."

        lines = generate_igv_batch_per_sample(
            regions, bam_list, output, max_height, snapshot_dir, genome_build,
        )

        all_lines.extend(lines)

    footer = generate_igv_batch_footer()
    all_lines.extend(footer)

    output_lines(all_lines, output)


def close_files(args):
    """Close all files in given iterable."""
    args_dict = vars(args)
    for arg_value in args_dict.values():
        if hasattr(arg_value, "close"):
            arg_value.close()


if __name__ == "__main__":
    """Run the main function only if this file is run as a script."""
    main()
