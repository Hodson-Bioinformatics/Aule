#!/usr/bin/env snakemake


##### ATTRIBUTION #####

# Original Author:  Ilias Moutsopoulos
# Module Author:    Ilias Moutsopoulos
# Contributors:     Natasha Cutmore


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["aule_out_tables"]`
CFG = op.setup_module(
    name = "aule_out_tables",
    version = "1.0",
    subdirectories = ["inputs", "by_patient", "by_sample", "by_mutation", "by_sv", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _aule_out_tables_input_tables,
    _aule_out_tables_input_vcf,
    _aule_out_tables_input_maf,
    # _aule_out_tables_input_sv_merge


##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _aule_out_tables_input_tables:
    input:
        samples = CFG["inputs"]["samples"],
        patient_info = CFG["inputs"]["patient_info"],
        sample_info = CFG["inputs"]["sample_info"],
        whitelist = CFG["inputs"]["whitelist"],
        drivers = CFG["inputs"]["drivers"]
    output:
        samples = CFG["dirs"]["inputs"] + "metadata/{seq_type}--{genome_build}/samples.tsv",
        patient_info = CFG["dirs"]["inputs"] + "metadata/{seq_type}--{genome_build}/patient_info.csv",
        sample_info = CFG["dirs"]["inputs"] + "metadata/{seq_type}--{genome_build}/sample_info.csv",
        whitelist = CFG["dirs"]["inputs"] + "metadata/{seq_type}--{genome_build}/whitelist.csv",
        drivers = CFG["dirs"]["inputs"] + "metadata/{seq_type}--{genome_build}/drivers.txt",
    run:
        op.absolute_symlink(input.samples, output.samples)
        op.absolute_symlink(input.patient_info, output.patient_info)
        op.absolute_symlink(input.sample_info, output.sample_info)
        op.absolute_symlink(input.whitelist, output.whitelist)
        op.absolute_symlink(input.drivers, output.drivers)

rule _aule_out_tables_input_vcf:
    input:
        vcf = CFG["inputs"]["sample_vcf"]
    output:
        vcf = CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.vcf.gz"
    run:
        op.absolute_symlink(input.vcf, output.vcf)

rule _aule_out_tables_input_maf:
    input:
        maf = CFG["inputs"]["sample_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)

# rule _aule_out_tables_input_sv_merge:
#     input:
#         csv = CFG["inputs"]["sv_csv"]
#     output:
#         csv = CFG["dirs"]["inputs"] + "sv_merge/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sv_merge_output.csv"
#     run:
#         op.absolute_symlink(input.csv, output.csv)

rule _aule_out_tables_by_patient:
    input: 
        samples = str(rules._aule_out_tables_input_tables.output.samples),
        patient_info = str(rules._aule_out_tables_input_tables.output.patient_info),
    output:
        csv = CFG["dirs"]["by_patient"] + "{seq_type}--{genome_build}/by_patient.csv"
    params:
        script = CFG["scripts"]["join_sample_info"]
    log:
        stdout = CFG["logs"]["by_patient"] + "{seq_type}--{genome_build}/by_patient.stdout.log",
        stderr = CFG["logs"]["by_patient"] + "{seq_type}--{genome_build}/by_patient.stderr.log"
    conda:
        CFG["conda_envs"]["aule_out_tables"]
    resources:
        **CFG["resources"]["aule_out_tables"]
    shell:
        op.as_one_line("""
        Rscript --vanilla
                {params.script}
                --samples {input.samples}
                --in_csv {input.patient_info}
                --out_csv {output.csv}
                --column_name patient_id
                > {log.stdout} 2> {log.stderr}
        """)

rule _aule_out_tables_by_sample:
    input: 
        samples = str(rules._aule_out_tables_input_tables.output.samples),
        sample_info = str(rules._aule_out_tables_input_tables.output.sample_info),
    output:
        csv = CFG["dirs"]["by_sample"] + "{seq_type}--{genome_build}/by_sample.csv"
    params:
        script = CFG["scripts"]["join_sample_info"]
    log:
        stdout = CFG["logs"]["by_sample"] + "{seq_type}--{genome_build}/by_sample.stdout.log",
        stderr = CFG["logs"]["by_sample"] + "{seq_type}--{genome_build}/by_sample.stderr.log"
    conda:
        CFG["conda_envs"]["aule_out_tables"]
    resources:
        **CFG["resources"]["aule_out_tables"]
    shell:
        op.as_one_line("""
        Rscript --vanilla
                {params.script}
                --samples {input.samples}
                --in_csv {input.sample_info}
                --out_csv {output.csv}
                --column_name CMDL_ID
                > {log.stdout} 2> {log.stderr}
        """)

rule _aule_out_tables_by_mutation:
    input:
        vcfs = expand(
                [
                    str(rules._aule_out_tables_input_vcf.output.vcf)
                ],
                zip,  # Run expand() with zip(), not product()
                seq_type=CFG["runs"]["tumour_seq_type"],
                genome_build=CFG["runs"]["tumour_genome_build"],
                tumour_id=CFG["runs"]["tumour_sample_id"],
                normal_id=CFG["runs"]["normal_sample_id"],
                pair_status=CFG["runs"]["pair_status"]
        ),
        mafs = expand(
            expand(
                [
                    str(rules._aule_out_tables_input_maf.output.maf)
                ],
                zip,  # Run expand() with zip(), not product()
                seq_type=CFG["runs"]["tumour_seq_type"],
                genome_build=CFG["runs"]["tumour_genome_build"],
                tumour_id=CFG["runs"]["tumour_sample_id"],
                normal_id=CFG["runs"]["normal_sample_id"],
                pair_status=CFG["runs"]["pair_status"],
                allow_missing = True
            ),
            filter = CFG["options"]["maf_filter"],
            maf = CFG["options"]["maf_type"],
            base_name = CFG["options"]["maf_base_name"]
        ),
        whitelist = str(rules._aule_out_tables_input_tables.output.whitelist),
        drivers = str(rules._aule_out_tables_input_tables.output.drivers),
        sample_info = str(rules._aule_out_tables_input_tables.output.sample_info)
    output:
        csv = CFG["dirs"]["by_mutation"] + "{seq_type}--{genome_build}/by_mutation.csv",
        csv_clean = CFG["dirs"]["by_mutation"] + "{seq_type}--{genome_build}/by_mutation_clean.csv"
    params:
        vaf_filter = CFG["options"]["vaf_filter"],
        shm_distance = CFG["options"]["shm_distance"]
    log:
        stdout = CFG["logs"]["by_mutation"] + "{seq_type}--{genome_build}/by_mutation.stdout.log",
        stderr = CFG["logs"]["by_mutation"] + "{seq_type}--{genome_build}/by_mutation.stderr.log"
    conda:
        CFG["conda_envs"]["aule_out_tables"]
    resources:
        **CFG["resources"]["aule_out_tables"]
    script:
        "etc/R/add_mutation_info.R"

# rule _aule_out_tables_by_sv:
#     input:
#         csv = expand(
#                 [
#                     str(rules._aule_out_tables_input_sv_merge.output.csv)
#                 ],
#                 zip,  # Run expand() with zip(), not product()
#                 seq_type=CFG["runs"]["tumour_seq_type"],
#                 genome_build=CFG["runs"]["tumour_genome_build"],
#                 tumour_id=CFG["runs"]["tumour_sample_id"],
#                 normal_id=CFG["runs"]["normal_sample_id"],
#                 pair_status=CFG["runs"]["pair_status"]
#         ),
#         sample_info = str(rules._aule_out_tables_input_tables.output.sample_info)
#     output:
#         csv_all = CFG["dirs"]["by_sv"] + "{seq_type}--{genome_build}/by_sv_all.csv",
#         csv_roi = CFG["dirs"]["by_sv"] + "{seq_type}--{genome_build}/by_sv_sample.csv",
#         csv_simple = CFG["dirs"]["by_sv"] + "{seq_type}--{genome_build}/by_sv_patient.csv"
#     params:
#         join_by = CFG["options"]["column_name"]
#     log:
#         stdout = CFG["logs"]["by_sv"] + "{seq_type}--{genome_build}/by_sv.stdout.log",
#         stderr = CFG["logs"]["by_sv"] + "{seq_type}--{genome_build}/by_sv.stderr.log"
#     conda:
#         CFG["conda_envs"]["aule_out_tables"]
#     script:
#         "etc/R/combine_sv2.R"


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _aule_out_tables_output_tables:
    input:
        by_patient = str(rules._aule_out_tables_by_patient.output.csv),
        by_sample = str(rules._aule_out_tables_by_sample.output.csv),
        by_mutation = str(rules._aule_out_tables_by_mutation.output.csv),
        by_mutation_clean = str(rules._aule_out_tables_by_mutation.output.csv_clean),
        # by_sv_all = str(rules._aule_out_tables_by_sv.output.csv_all),
        # by_sv_sample = str(rules._aule_out_tables_by_sv.output.csv_roi),
        # by_sv_patient = str(rules._aule_out_tables_by_sv.output.csv_simple)
    output:
        by_patient = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/by_patient.csv",
        by_sample = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/by_sample.csv",
        by_mutation = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/by_mutation.csv",
        by_mutation_clean = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/by_mutation_clean.csv",
        # by_sv_all = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/by_sv_all.csv",
        # by_sv_sample = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/by_sv_sample.csv",
        # by_sv_patient = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/by_sv_patient.csv"
    run:
        op.relative_symlink(input.by_patient, output.by_patient, in_module = True)
        op.relative_symlink(input.by_sample, output.by_sample, in_module = True)
        op.relative_symlink(input.by_mutation, output.by_mutation, in_module = True)
        op.relative_symlink(input.by_mutation_clean, output.by_mutation_clean, in_module = True)
        # op.relative_symlink(input.by_sv_all, output.by_sv_all, in_module = True)
        # op.relative_symlink(input.by_sv_sample, output.by_sv_sample, in_module = True)
        # op.relative_symlink(input.by_sv_patient, output.by_sv_patient, in_module = True)
 
 
# # Generates the target sentinels for each run, which generate the symlinks
rule _aule_out_tables_all:
    input:
        expand(
            [
                str(rules._aule_out_tables_output_tables.output.by_patient),
                str(rules._aule_out_tables_output_tables.output.by_sample),
                str(rules._aule_out_tables_output_tables.output.by_mutation),
                str(rules._aule_out_tables_output_tables.output.by_mutation_clean),
                # str(rules._aule_out_tables_output_tables.output.by_sv_all),
                # str(rules._aule_out_tables_output_tables.output.by_sv_sample),
                # str(rules._aule_out_tables_output_tables.output.by_sv_patient)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"]
        )

##### CLEANUP #####

# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
