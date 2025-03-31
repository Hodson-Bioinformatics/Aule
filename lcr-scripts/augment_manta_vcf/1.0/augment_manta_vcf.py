#!/usr/bin/env python3

"""
augment_manta_vcf.py
====================

This script augments Manta VCF files with additional information. The following
fields are added to the INFO and FORMAT/SAMPLE columns where applicable:

    1) TR (FORMAT field): Sum of the SR and PR FORMAT fields per allele.
    2) DP (FORMAT field): Sum of the SR and PR FORMAT fields across all alleles.
    3) VAF (FORMAT field): Variant allele fraction for the alternate allele.
    4) REGIONS (INFO field): List of regions from the given BED files that
       overlap each variant position. One example use case is labelling
       breakpoints occuring in loci of interest such as MYC, BCL2, BCL6,
       and the three IG loci. This can readily be done downstream, but
       having the breakpoints labelled in the VCF file can enable quick
       exploratory analyses.

The script can optionally update the sample IDs in the VCF header, namely
the row with "#CHROM". Unfortunately, different Manta output VCF files
include different combinations of the tumour and/or normal samples used
in the analysis. Based on existing Manta output files, it was found that
the sample present in the VCF files follow the rules below. The script
will attempt to infer the "VCF type" (e.g., "somaticSV") based on the
input VCF file path. If this isn't possible for whatever reason, the
user can provide a value for the `--vcf_type` argument. Depending on
the VCF type, the user will need to provide values for `--tumour_id`
and/or `--normal_id`.

    1) candidateSV: No samples
    2) candidateSmallIndels: No samples
    3) rnaSV: Tumour sample only
    4) tumorSV: Tumour sample only
    5) diploidSV: Normal sample only
    6) somaticSV: Normal sample first and then tumour sample

Inputs
------
- Manta VCF file

Outputs
-------
- Augmented VCF file

Caveats
-------
- This script assumes that the Manta VCF files only contain one alternate
  allele, even though the VCF header hints that there could be more than
  one (e.g., "for the ref and alt alleles in the order listed"). Out of
  hundreds of Manta VCF files, none had more than one alternate allele,
  so support for more alternate allele won't be added for now.
- Due to limitations with the `cyvcf2` API, it is not possible to update
  the sample IDs during the pass when the variants are being processed.
  Unfortunately, a second pass of the VCF file is required. Updating the
  sample IDs is done as the second step to avoid dealing with compressed
  files since cyvcf2 can more easily handle gzip-compressed VCF files.
"""


# Import standard modules
import os
import sys
import shutil
import warnings
import argparse
import tempfile

# Import third-party modules
import pyranges as pr
import numpy as np
from cyvcf2 import VCF, Writer


# Track verbosity globally (to avoid countless function parameters)
IS_QUIET = False


def main():
    """Runs all of the other functions."""

    # Parse command-line arguments
    args = parse_arguments()

    # Change verbosity
    global IS_QUIET
    IS_QUIET = args.quiet

    # Add new fields to each variant
    augment_vcf(args.vcf_input, args.vcf_output, args.bed_regions, args.decimals)

    # Parse VCF type if not provided
    if args.vcf_type is None:
        args.vcf_type = parse_vcf_type(args.vcf_input)

    # Update sample IDs (will skip if both are None)
    update_sample_ids(args.vcf_output, args.vcf_type, args.tumour_id, args.normal_id)


def parse_arguments():
    """Parses and validates command-line arguments."""

    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf_input", help="Manta VCF file (gzip-compressed or not).")
    parser.add_argument("vcf_output", help="Output (augmented) VCF file.")
    parser.add_argument(
        "--bed_regions",
        "-b",
        nargs="+",
        default=[],
        help="Regions BED file(s). Names (column 4) must be unique across all files.",
    )
    parser.add_argument(
        "--decimals",
        "-d",
        type=int,
        default=3,
        help="Number of decimals for rounding.",
    )
    parser.add_argument("--tumour_id", "-t", help="Tumour sample ID.")
    parser.add_argument("--normal_id", "-n", help="Normal sample ID.")
    parser.add_argument(
        "--vcf_type",
        "-v",
        choices=[
            "tumorSV",
            "rnaSV",
            "somaticSV",
            "diploidSV",
            "candidateSmallIndels",
            "candidateSV",
        ],
        help="Manta VCF type.",
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="Silence warnings.",
    )
    args = parser.parse_args()

    # Validate command-line arguments
    assert not args.vcf_output.endswith(
        ".gz"
    ), "Cannot output gzip-compressed VCF file."

    return args


def modify_header(vcf, bed_files):
    """Adds info on TR and VAF FORMAT fields to VCF header."""

    # Add FORMAT and INFO field information to header
    vcf.add_format_to_header(
        {
            "ID": "TR",
            "Description": "Total read support (SR + PR) for the ref and alt alleles in the order listed",
            "Type": "Float",
            "Number": ".",
        }
    )
    vcf.add_format_to_header(
        {
            "ID": "DP",
            "Description": "Total read support (SR + PR) for the ref and first alt allele combined",
            "Type": "Float",
            "Number": "1",
        }
    )
    vcf.add_format_to_header(
        {
            "ID": "VAF",
            "Description": "Variant allele fraction for the first alt allele",
            "Type": "Float",
            "Number": "1",
        }
    )

    # Only add info about REGIONS field if at least one BED is given
    if len(bed_files) > 0:
        vcf.add_info_to_header(
            {
                "ID": "REGIONS",
                "Description": (
                    "Overlapping regions from the given BED files "
                    "(see `##regions_bed` header lines)"
                ),
                "Type": "String",
                "Number": ".",
            }
        )

    # Add list of regions BED files to header
    for bed_file in bed_files:
        bed_file = os.path.abspath(bed_file)
        vcf.add_to_header(f"##regions_bed={bed_file}\n")

    # Add command to header
    sys.argv[0] = os.path.abspath(sys.argv[0])
    command = " ".join(sys.argv)
    vcf.add_to_header(f"##cmdline={command}\n")

    return vcf


def parse_bed_files(bed_files):
    """Creates PyRanges objects from the BED files."""

    # Skip if no BED files are provided
    if len(bed_files) == 0:
        return

    # Load BED files
    beds = [pr.read_bed(b) for b in bed_files]

    # Check that all BED files have the first four columns
    for bed_file, bed in zip(bed_files, beds):
        assert "Name" in bed.columns, f"Name (column 4) missing from {bed_file}."

    # Concatenate BED files and only keep Name column
    bed = pr.concat(beds)
    bed = bed.unstrand()
    bed = bed[["Name"]]

    # Ensure unique names
    assert bed.Name.is_unique, "Names (column 4) not unique across BED files."

    return bed


def get_overlapping_region_names(variant, bed):
    """Gets list of overlapping region names from a BED file."""

    # Generate 0-based region for querying regions
    sv_type = variant.INFO["SVTYPE"]
    chrom = variant.CHROM
    start = variant.POS if sv_type == "INS" else variant.POS - 1
    end = variant.POS if sv_type == "BND" else variant.INFO["END"]
    end = end - 1 if sv_type == "INS" else end

    # Get overlapping regions
    overlapping_regions = bed[chrom, start:end]

    # Return empty list if no overlapping regions are found
    if len(overlapping_regions) == 0:
        return []

    # Issue warning if there are many overlapping regions
    if len(overlapping_regions) > 10 and not IS_QUIET:
        warnings.warn(
            "Warning: More than 10 overlapping regions. If this is "
            "expected, you can silence this warning with `--quiet`."
        )

    return overlapping_regions.Name.tolist()


def add_fields_to_variant(variant, bed, decimals):
    """Adds the given variant with additional fields.

    The new fields are described in the top-level docstring.
    """

    # Extract split read (sr) and read pair (pr) counts
    sr = variant.format("SR")
    pr = variant.format("PR")

    # Fill in any missing fields with zeros
    if sr is None or pr is None:
        # Use the same shape as the non-missing field
        shape = sr.shape if pr is None else pr.shape
        zeros = np.zeros(shape)
        sr = zeros if sr is None else sr
        pr = zeros if pr is None else pr

    # Calculate total read counts (tr)
    tr = sr + pr

    # Expect two alleles (one REF and one ALT)
    num_alleles = tr.shape[1]
    assert num_alleles == 2, "Encountered more than one alternate allele."

    # Calculate allele-specific counts, total depth (dp), and VAF
    ref_count = tr[:, 0]
    alt_count = tr[:, 1]
    dp = ref_count + alt_count
    vaf = alt_count / dp
    vaf = vaf.round(decimals)

    # Set the calculated values as FORMAT fields
    variant.set_format("TR", tr)
    variant.set_format("DP", dp)
    variant.set_format("VAF", vaf)

    # Get list of overlapping region names and add if not empty
    if bed is not None:
        region_names = get_overlapping_region_names(variant, bed)
        if len(region_names) > 0:
            variant.INFO["REGIONS"] = ",".join(region_names)

    return variant


def augment_vcf(vcf_in_file, vcf_out_file, bed_files, decimals):
    """Parses and augments VCF file."""

    # Read in the input VCF file
    vcf_in = VCF(vcf_in_file)

    # Add rows to the header for each new field
    vcf_in = modify_header(vcf_in, bed_files)

    # Set up a write based on the tweaked input VCF file
    vcf_out = Writer(vcf_out_file, vcf_in)

    # Parse BED files
    bed = parse_bed_files(bed_files)

    # Iterate over every variant record
    for variant in vcf_in:
        # Augment the variant by adding new fields (if there are samples)
        num_samples = len(vcf_in.samples)
        if num_samples > 0:
            variant = add_fields_to_variant(variant, bed, decimals)
        # Output the augmented variant
        vcf_out.write_record(variant)

    # Close input and output VCF files
    vcf_in.close()
    vcf_out.close()


def parse_vcf_type(filename):
    """Extracts VCF type from the file path."""

    # Create list of all possible Manta VCF types
    vcf_types = [
        "candidateSV",
        "candidateSmallIndels",
        "diploidSV",
        "somaticSV",
        "tumorSV",
        "rnaSV",
    ]

    # Filter down for those that appear in the filename
    vcf_types = filter(lambda x: x in filename, vcf_types)
    vcf_types = list(vcf_types)
    num_matches = len(vcf_types)

    # Check that only one match exists
    assert num_matches == 1, (
        "Could not infer VCF type from file name. Expected 1 match, but found "
        f"{num_matches} instead. Please use `--vcf_type` to disambiguate."
    )

    return vcf_types[0]


def update_header_line(line, vcf_type, tumour_id, normal_id):
    """Updates header line based on VCF type.

    See top-level docstring for more details.
    """

    # Split line into columns
    columns = line.rstrip("\n").split("\t")
    num_columns = len(columns)

    # No samples
    if vcf_type in ["candidateSV", "candidateSmallIndels"]:
        assert num_columns == 8, f"Expected 8 columns, found {num_columns}."

    # Tumour sample only
    elif vcf_type in ["rnaSV", "tumorSV"]:
        assert num_columns == 10, f"Expected 10 columns, found {num_columns}."
        assert tumour_id is not None, f"`--tumour_id` is required for {vcf_type}."
        columns[9] = tumour_id

    # Normal sample only
    elif vcf_type in ["diploidSV"]:
        assert num_columns == 10, f"Expected 10 columns, found {num_columns}."
        assert normal_id is not None, f"`--normal_id` is required for {vcf_type}."
        columns[9] = normal_id

    # Normal sample first and then tumour sample
    elif vcf_type in ["somaticSV"]:
        assert num_columns == 11, f"Expected 11 columns, found {num_columns}."
        assert tumour_id is not None, f"`--tumour_id` is required for {vcf_type}."
        assert normal_id is not None, f"`--normal_id` is required for {vcf_type}."
        columns[9] = normal_id
        columns[10] = tumour_id

    # This shouldn't happen because either it was provided by the user, which
    # is restricted to the choices given to `argparse`, or it was parsed by
    # `parse_vcf_type`, which shouldn't return anything unexpected.
    else:
        raise Exception(f"Encountered unknown VCF type ({vcf_type}).")

    return "\t".join(columns) + "\n"


def update_sample_ids(vcf_file, vcf_type, tumour_id, normal_id):
    """Processes VCF file to update sample IDs in header."""

    # Skip if both tumour_id and normal_id are None
    if tumour_id is None and normal_id is None:
        return

    # Open files
    _, temp_file = tempfile.mkstemp()
    with open(vcf_file, "r") as vcf, open(temp_file, "w") as temp:

        # Iterate over every line (and only process header)
        for line in vcf:
            # Find the header line
            if line.startswith("#CHROM"):
                line = update_header_line(line, vcf_type, tumour_id, normal_id)
            temp.write(line)

    # Replace the old VCF file with the updated file
    shutil.move(temp_file, vcf_file)


if __name__ == "__main__":
    main()
