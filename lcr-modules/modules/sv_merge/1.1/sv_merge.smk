#!/usr/bin/env snakemake

##### ATTRIBUTION #####

# Original Author:  Natasha Cutmore
# Module Author:    Ilias Moutsopoulos
# Contributors:     Natasha Cutmore

##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["sv_merge"]`
CFG = op.setup_module(
    name = "sv_merge",
    version = "1.1",
    subdirectories = ["inputs", "merge", "filter", "summarise", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _sv_merge_input_vcfs,
    _sv_merge_input_metadata,
    _sv_merge_output_tables,
    _sv_merge_all


##### RULES #####
# Symlinks the input vcf files into the module results directory (under '00-inputs/')
rule _sv_merge_input_vcfs:
    input:
        gridss_vcf = CFG["inputs"]["gridss_vcf"],
        manta_vcf = CFG["inputs"]["manta_vcf"]
    output:
        gridss_vcf = CFG["dirs"]["inputs"] + "GRIDSS/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.gridss_somatic_filtered.vcf.gz",
        manta_vcf = CFG["dirs"]["inputs"] + "Manta/{seq_type}--{genome_build}/somaticSV/{tumour_id}--{normal_id}--{pair_status}.somaticSV.vcf"
    run:
        op.absolute_symlink(input.gridss_vcf, output.gridss_vcf)
        op.absolute_symlink(input.manta_vcf, output.manta_vcf)

# Symlinks the input meta data files into the module results directory (under '00-inputs/')
rule _sv_merge_input_metadata:
    input:
        samples = CFG["inputs"]["samples"],
        roi_annotations = CFG["inputs"]["roi_annotations"]
    output:
        samples = CFG["dirs"]["inputs"] + "metadata/samples.tsv",
        roi_annotations = CFG["dirs"]["inputs"] + "metadata/roi_annotations.csv"
    run:
        op.absolute_symlink(input.samples, output.samples)
        op.absolute_symlink(input.roi_annotations, output.roi_annotations)

# Merge output from callers
rule _sv_merge_merge_caller_output:
    input:
        gridss_vcf = str(rules._sv_merge_input_vcfs.output.gridss_vcf),
        manta_vcf = str(rules._sv_merge_input_vcfs.output.manta_vcf)
    output:
        csv = CFG["dirs"]["merge"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_sv_merge_all_variants.csv"
    log:
        stdout = CFG["logs"]["merge"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.outputs.stdout.log",
        stderr = CFG["logs"]["merge"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.outputs.stderr.log"
    conda:
        CFG["conda_envs"]["postprocessing"]
    resources:
        **CFG["resources"]["_sv_merge_run"]
    script:
        "src/R/merge_caller_output.R"

# Filter merged output by regions of interest
rule _sv_merge_filter_output:
    input:
        merged_csv = str(rules._sv_merge_merge_caller_output.output.csv),
        roi_annotations = str(rules._sv_merge_input_metadata.output.roi_annotations)
    output:
        csv = CFG["dirs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_sv_merge_roi_filtered.csv"
    log:
        stdout = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.outputs.stdout.log",
        stderr = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.outputs.stderr.log"
    conda:
        CFG["conda_envs"]["postprocessing"]
    resources:
        **CFG["resources"]["_sv_merge_run"]
    script:
        "src/R/roi_filter_output.R"

# Create summary table
rule _sv_merge_summarise_output:
    input:
        csvs = expand(
                [
                    str(rules._sv_merge_filter_output.output.csv)
                ],
                zip,  # Run expand() with zip(), not product()
                seq_type=CFG["runs"]["tumour_seq_type"],
                genome_build=CFG["runs"]["tumour_genome_build"],
                tumour_id=CFG["runs"]["tumour_sample_id"],
                normal_id=CFG["runs"]["normal_sample_id"],
                pair_status=CFG["runs"]["pair_status"],
                allow_missing = True
            ),
        samples = str(rules._sv_merge_input_metadata.output.samples)
    output:
        csv = CFG["dirs"]["summarise"] + "{seq_type}--{genome_build}/sv_merge_summary.csv"
    params:
        SR_filter = CFG["params"]["SR_filter"],
        RP_filter = CFG["params"]["RP_filter"]
    log:
        stdout = CFG["logs"]["summarise"] + "{seq_type}--{genome_build}/svs_summarised.outputs.stdout.log",
        stderr = CFG["logs"]["summarise"] + "{seq_type}--{genome_build}/svs_summarised.outputs.stderr.log"
    conda:
        CFG["conda_envs"]["postprocessing"]
    resources:
        **CFG["resources"]["_sv_merge_run"]
    script:
        "src/R/summarise_output.R"

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _sv_merge_output_tables:
    input:
        csv = str(rules._sv_merge_summarise_output.output.csv),
    output:
        csv = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/sv_merge_summary.csv",
    run:
        op.relative_symlink(input.csv, output.csv, in_module= True)

rule _sv_merge_all:
    input:
        expand(
            [
                str(rules._sv_merge_output_tables.output.csv)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"]
            )

##### CLEANUP #####

# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
