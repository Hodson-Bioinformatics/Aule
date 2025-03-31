#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ilias Moutsopoulos
# Module Author:    Ilias Moutsopoulos
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op

# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version="1.0.11"
import pkg_resources
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

# To avoid this we need to add the "packaging" module as a dependency for LCR-modules or oncopipe

current_version = pkg_resources.get_distribution("oncopipe").version
if version.parse(current_version) < version.parse(min_oncopipe_version):
    print('\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}' + '\x1b[0m')
    print('\x1b[0;31;40m' + f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m')
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["aule_error"]`
CFG = op.setup_module(
    name = "aule_error",
    version = "1.0",
    subdirectories = ["inputs", "bed_to_coords", "collapse_errors", "merge_sample_errors", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _aule_error_input_bam,
    _aule_error_input_bed,
    _aule_error_output_csv,
    _aule_error_all,


##### RULES #####


# Symlink the input files into the module results directory (under '00-inputs/')
rule _aule_error_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        op.absolute_symlink(input.bam, output.bam)


rule _aule_error_input_bed:
    input:
        bed = CFG["inputs"]["regions_bed"]
    output:
        bed = CFG["dirs"]["inputs"] + "bed/{seq_type}--{genome_build}/regions.bed"
    run:
        op.absolute_symlink(input.bed, output.bed)


rule _aule_error_bed_to_coords:
    input:
        bed = str(rules._aule_error_input_bed.output.bed)
    output:
        coords = CFG["dirs"]["bed_to_coords"] + "{seq_type}--{genome_build}/coords.csv"
    log:
        stdout = CFG["logs"]["bed_to_coords"] + "{seq_type}--{genome_build}/bed_to_coords.stdout.log",
        stderr = CFG["logs"]["bed_to_coords"] + "{seq_type}--{genome_build}/bed_to_coords.stderr.log"
    params:
        script = CFG["scripts"]["bed_to_coords"]
    conda:
        CFG["conda_envs"]["bed_to_coords"]
    shell:
        op.as_one_line("""
        RScript --vanilla
                {params.script}
                --in_bed {input.bed}
                --out_csv {output.coords}
                > {log.stdout} 2> {log.stderr}
        """)


rule _aule_error_collapse_errors:
    input:
        bam = str(rules._aule_error_input_bam.output.bam),
        coords = str(rules._aule_error_bed_to_coords.output.coords)
        # need to confirm bed2coords conda env works
    output:
        errors = CFG["dirs"]["collapse_errors"] + "{seq_type}--{genome_build}/{sample_id}/errors_collapsed.csv",
        stats = CFG["dirs"]["collapse_errors"] + "{seq_type}--{genome_build}/{sample_id}/stats.csv"
    log:
        stdout = CFG["logs"]["collapse_errors"] + "{seq_type}--{genome_build}/{sample_id}/collapse_errors.stdout.log",
        stderr = CFG["logs"]["collapse_errors"] + "{seq_type}--{genome_build}/{sample_id}/collapse_errors.stderr.log"
    params:
        script = CFG["scripts"]["collapse_errors"]
    conda:
        CFG["conda_envs"]["collapse_errors"]
    shell:
        op.as_one_line("""
        python {params.script}
               -in_bam {input.bam}
               -in_coords {input.coords}
               -out_errors {output.errors}
               -out_stats {output.stats}
               > {log.stdout} 2> {log.stderr}
        """)


rule _aule_error_merge_sample_errors:
    input:
        errors = expand(
                [
                    str(rules._aule_error_collapse_errors.output.errors)
                ],
                zip,  # Run expand() with zip(), not product()
                seq_type=CFG["samples"]["seq_type"],
                genome_build=CFG["samples"]["genome_build"],
                sample_id=CFG["samples"]["sample_id"]
        ),
        stats = expand(
                [
                    str(rules._aule_error_collapse_errors.output.stats)
                ],
                zip,  # Run expand() with zip(), not product()
                seq_type=CFG["samples"]["seq_type"],
                genome_build=CFG["samples"]["genome_build"],
                sample_id=CFG["samples"]["sample_id"]
        )
    output:
        errors_split_by_sample = CFG["dirs"]["merge_sample_errors"] + "{seq_type}--{genome_build}/errors_split_by_sample.csv",
        errors_merged = CFG["dirs"]["merge_sample_errors"] + "{seq_type}--{genome_build}/errors_merged.csv",
        stats = CFG["dirs"]["merge_sample_errors"] + "{seq_type}--{genome_build}/stats.csv"
    log:
        stdout = CFG["logs"]["merge_sample_errors"] + "{seq_type}--{genome_build}/merge_sample_errors.stdout.log",
        stderr = CFG["logs"]["merge_sample_errors"] + "{seq_type}--{genome_build}/merge_sample_errors.stderr.log"
    params:
        script = CFG["scripts"]["merge_sample_errors"]
    conda:
        CFG["conda_envs"]["merge_sample_errors"]
    script:
        "{params.script}" # "etc/R/merge_sample_errors.R"


rule _aule_error_output_csv:
    input:
        errors_split_by_sample = str(rules._aule_error_merge_sample_errors.output.errors_split_by_sample),
        errors_merged = str(rules._aule_error_merge_sample_errors.output.errors_merged),
        stats = str(rules._aule_error_merge_sample_errors.output.stats)
    output:
        errors_split_by_sample = CFG["dirs"]["outputs"] + "csv/{seq_type}--{genome_build}/errors_split_by_sample.csv",
        errors_merged = CFG["dirs"]["outputs"] + "csv/{seq_type}--{genome_build}/errors_merged.csv",
        stats = CFG["dirs"]["outputs"] + "csv/{seq_type}--{genome_build}/stats.csv"
    run:
        op.relative_symlink(input.errors_split_by_sample, output.errors_split_by_sample, in_module = True)
        op.relative_symlink(input.errors_merged, output.errors_merged, in_module = True)
        op.relative_symlink(input.stats, output.stats, in_module = True)


# Generates the target sentinels for each run, which generate the symlinks
rule _aule_error_all:
    input:
        expand(
            [
                str(rules._aule_error_output_csv.output.errors_split_by_sample),
                str(rules._aule_error_output_csv.output.errors_merged),
                str(rules._aule_error_output_csv.output.stats),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
