#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Natasha Cutmore
# Module Author:    Natasha Cutmore
# Contributors:     N/A

# This snakemake file produces output that is formatted correctly for input #into the NCI online LymphGen predictor tool.  

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
# `CFG` is a shortcut to `config["lcr-modules"]["lymphgen_subtyping"]`
CFG = op.setup_module(
    name = "lymphgen_subtyping",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "preproc_input", "lymphgen_subtyping", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _lymphgen_subtyping_input,
    _lymphgen_subtyping_sample_annotation_file
    _lymphgen_subtyping_output_csv,
    _lymphgen_subtyping_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
# TODO: If applicable, add an input rule for each input file used by the module
# TODO: If applicable, create second symlink to .crai file in the input function, to accomplish cram support
rule _lymphgen_subtyping_input:
    input:
        by_sample = CFG["inputs"]["by_sample"],
        by_mutation = CFG["inputs"]["by_mutation"],
        by_sv_samples = CFG["inputs"]["by_sv_samples"]
    output:
        by_sample = CFG["dirs"]["inputs"] + "csv/{seq_type}--{genome_build}/by_sample.csv",
        by_mutation = CFG["dirs"]["inputs"] + "csv/{seq_type}--{genome_build}/by_mutation.csv",
        by_sv_samples = CFG["dirs"]["inputs"] + "csv/{seq_type}--{genome_build}/by_sv_samples.csv"
    run:
        op.absolute_symlink(input.by_sample, output.by_sample)
        op.absolute_symlink(input.by_mutation, output.by_mutation)
        op.absolute_symlink(input.by_sv_samples, output.by_sv_samples)


rule _lymphgen_subtyping_sample_annotation_file:
    input:
        by_sv_samples = str(rules._lymphgen_subtyping_input.output.by_sv_samples),
        by_sample = str(rules._lymphgen_subtyping_input.output.by_sample)
    output:
        sample_annotation_file_txt = CFG["dirs"]["preproc_input"] + "{seq_type}--{genome_build}/sample_annotation_file.txt"
    params:
        script = CFG["scripts"]["sample_annotation_file"]
    log:
        stdout = CFG["logs"]["preproc_input"] + "{seq_type}--{genome_build}/sample_annotation_file.stdout.log",
        stderr = CFG["logs"]["preproc_input"] + "{seq_type}--{genome_build}/sample_annotation_file.stderr.log"
    conda:
        CFG["conda_envs"]["postprocessing"]
    shell:
        op.as_one_line("""
        Rscript --vanilla
                {params.script}
                --by_sv_samples {input.by_sv_samples}
                --samples {input.by_sample}
                --out_txt {output.sample_annotation_file_txt}
                > {log.stdout} 2> {log.stderr}
        """)

rule _lymphgen_subtyping_mutation_flat_file:
    input:
        by_sample = str(rules._lymphgen_subtyping_input.output.by_sample),
        by_mutation = str(rules._lymphgen_subtyping_input.output.by_mutation)
    output:
        mutation_flat_file_txt = CFG["dirs"]["preproc_input"] + "{seq_type}--{genome_build}/mutation_flat_file.txt",
        mutation_gene_list_txt = CFG["dirs"]["preproc_input"] + "{seq_type}--{genome_build}/mutation_gene_list.txt"
    params:
        script = CFG["scripts"]["mutation_flat_file"]
    log:
        stdout = CFG["logs"]["preproc_input"] + "{seq_type}--{genome_build}/mutation_files.stdout.log",
        stderr = CFG["logs"]["preproc_input"] + "{seq_type}--{genome_build}/mutation_files.stderr.log"
    conda:
        CFG["conda_envs"]["postprocessing"]
    shell:
        op.as_one_line("""
        Rscript --vanilla
                {params.script}
                --samples {input.by_sample}
                --mutation_file {input.by_mutation}
                --out_mff {output.mutation_flat_file_txt}
                --out_mgl {output.mutation_gene_list_txt}
                > {log.stdout} 2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _lymphgen_subtyping_output:
    input:
        sample_annotation_file_txt = str(rules._lymphgen_subtyping_sample_annotation_file.output.sample_annotation_file_txt),
        mutation_flat_file_txt = str(rules._lymphgen_subtyping_mutation_flat_file.output.mutation_flat_file_txt),
        mutation_gene_list_txt = str(rules._lymphgen_subtyping_mutation_flat_file.output.mutation_gene_list_txt)
    output:
        sample_annotation_file_txt = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/sample_annotation_file.txt",
        mutation_flat_file_txt = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/mutation_flat_file.txt",
        mutation_gene_list_txt = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/mutation_gene_list.txt"
    run:
        op.relative_symlink(input.sample_annotation_file_txt, output.sample_annotation_file_txt, in_module= True),
        op.relative_symlink(input.mutation_flat_file_txt, output.mutation_flat_file_txt, in_module= True),
        op.relative_symlink(input.mutation_gene_list_txt, output.mutation_gene_list_txt, in_module= True)



# Generates the target sentinels for each run, which generate the symlinks
rule _lymphgen_subtyping_all:
    input:
        expand(
            [
                str(rules._lymphgen_subtyping_output.output.sample_annotation_file_txt),
                str(rules._lymphgen_subtyping_output.output.mutation_flat_file_txt),
                str(rules._lymphgen_subtyping_output.output.mutation_gene_list_txt)

                # TODO: If applicable, add other output rules here
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
