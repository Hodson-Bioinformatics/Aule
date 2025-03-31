#!/usr/bin/env snakemake


##### ATTRIBUTION #####

# Original Author:  Ilias Moutsopoulos
# Module Author:    Joanna Krupka
# Contributors:     Ilias Moutsopoulos, Natasha Cutmore


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["aule_out_tables"]`
CFG = op.setup_module(
    name = "aule_out_tables",
    version = "1.2",
    subdirectories = ["inputs", "by_mutation", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _aule_out_tables_input_tables,
    _aule_out_tables_input_vcf,
    _aule_out_tables_input_maf,
    _aule_out_tables_output_tables,
    _aule_out_tables_all

##### RULES #####

# Symlinks the input metadata files into the module results directory (under '00-inputs/')
rule _aule_out_tables_input_tables:
    input:
        whitelist = CFG["inputs"]["whitelist"],
        drivers = CFG["inputs"]["drivers"]
    output:
        whitelist = CFG["dirs"]["inputs"] + "metadata/{seq_type}--{genome_build}/whitelist.csv",
        drivers = CFG["dirs"]["inputs"] + "metadata/{seq_type}--{genome_build}/drivers.txt",
    run:
        op.absolute_symlink(input.whitelist, output.whitelist)
        op.absolute_symlink(input.drivers, output.drivers)

rule _aule_out_tables_input_blacklisted_coord:
    input:
        blacklisted_coord = CFG["inputs"]["blacklisted_coord"]
    output:
        blacklisted_coord = CFG["dirs"]["inputs"] + "blacklist/{seq_type}--{genome_build}/blacklisted_coord.csv",
    run:
        op.absolute_symlink(input.blacklisted_coord, output.blacklisted_coord)
        
rule _aule_out_tables_input_panel:
    input:
        panel = CFG["inputs"]["panel"]
    output:
        panel = CFG["dirs"]["inputs"] + "panel/{seq_type}--{genome_build}/panel.bed",
    run:
        op.absolute_symlink(input.panel, output.panel)
        
# Symlinks the input maf files into the module results directory (under '00-inputs/')
rule _aule_out_tables_input_maf:
    input:
        maf = CFG["inputs"]["sample_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)

# Generates table with comprehensively annotated mutation data
rule _aule_out_tables_by_mutation:
    input:
        mafs = expand(
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
        whitelist = str(rules._aule_out_tables_input_tables.output.whitelist),
        drivers = str(rules._aule_out_tables_input_tables.output.drivers),
        panel = str(rules._aule_out_tables_input_panel.output.panel),
        blacklisted_coord = str(rules._aule_out_tables_input_blacklisted_coord.output.blacklisted_coord)
    output:
        csv = CFG["dirs"]["by_mutation"] + "{seq_type}--{genome_build}/by_mutation.csv",
    params:
        shm_distance = CFG["options"]["shm_distance"],
        min_pon = CFG["options"]["min_pon"],
        filtering = CFG["options"]["filtering"],
        min_germline = CFG["options"]["min_germline"],
        max_error_prob = CFG["options"]["max_error_prob"],    
        min_umi = CFG["options"]["min_umi"]     
    log:
        stdout = CFG["logs"]["by_mutation"] + "{seq_type}--{genome_build}/by_mutation.stdout.log",
        stderr = CFG["logs"]["by_mutation"] + "{seq_type}--{genome_build}/by_mutation.stderr.log"
    conda:
        CFG["conda_envs"]["postprocessing"]
    threads:
        CFG["threads"]["by_mutation"]
    resources:
        **CFG["resources"]["by_mutation"]
    script:
        "etc/R/add_mutation_info.R"

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _aule_out_tables_output_tables:
    input:
        by_mutation = str(rules._aule_out_tables_by_mutation.output.csv),
    output:
        by_mutation = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/by_mutation.csv",
    run:
        op.relative_symlink(input.by_mutation, output.by_mutation, in_module = True)

# Generates the target sentinels for each run, which generate the symlinks
rule _aule_out_tables_all:
    input:
        expand(
            [
                str(rules._aule_out_tables_output_tables.output.by_mutation),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"]
        )

##### CLEANUP #####

# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)

