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
    version = "1.0",
    subdirectories = ["inputs", "filter_gridss", "filter_manta", "survivor", "merge_format", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _sv_merge_filter_gridss,
    _sv_merge_filter_manta,
    _sv_merge_prepare_survivor,
    _sv_merge_output_tables,
    _sv_merge_all


##### RULES #####
# Symlinks the input files into the module results directory (under '00-inputs/')
rule _sv_merge_inputs:
    input:
        igcaller_tsv = CFG["inputs"]["igcaller_tsv"],
        gridss_vcf = CFG["inputs"]["gridss_vcf"],
        manta_vcf = CFG["inputs"]["manta_vcf"]
    output:
        igcaller_tsv = CFG["dirs"]["inputs"] + "IgCaller/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.output.filtered.tsv",
        gridss_vcf = CFG["dirs"]["inputs"] + "GRIDSS/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.gridss_somatic_filtered.vcf.gz",
        manta_vcf = CFG["dirs"]["inputs"] + "Manta/{seq_type}--{genome_build}/somaticSV/{tumour_id}--{normal_id}--{pair_status}.somaticSV.vcf"
    run:
        op.absolute_symlink(input.igcaller_tsv, output.igcaller_tsv)
        op.absolute_symlink(input.gridss_vcf, output.gridss_vcf)
        op.absolute_symlink(input.manta_vcf, output.manta_vcf)


rule _sv_merge_filter_gridss:
    input: 
        vcf = str(rules._sv_merge_inputs.output.gridss_vcf),
        roi_coordinates = CFG["inputs"]["roi_coordinates"]
    output: 
        vcf_roi_filtered = CFG["dirs"]["filter_gridss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.gridss_roi_filtered.vcf"
    params:
        vaf_threshold = CFG["params"]["vaf_threshold"]
    shell:
        op.as_one_line("""
        bcftools view -T {input.roi_coordinates}
                      -i 'FILTER != "PON" & INFO/TAF>{params.vaf_threshold}' 
                      -O v 
                      -o {output.vcf_roi_filtered} 
                      {input.vcf}
        """)

rule _sv_merge_filter_manta:
    input: 
        vcf = str(rules._sv_merge_inputs.output.manta_vcf),
        roi_coordinates = CFG["inputs"]["roi_coordinates"]
    output: 
        vcf_roi_filtered = CFG["dirs"]["filter_manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.manta_roi_filtered.vcf"
    params:
        vaf_threshold = CFG["params"]["vaf_threshold"]
    shell:
        op.as_one_line("""
        bcftools view -T {input.roi_coordinates} 
                      -i 'FILTER = "PASS" & INFO/SOMATICSCORE>=50 & FMT/VAF>{params.vaf_threshold}' 
                      -O v 
                      -o {output.vcf_roi_filtered} 
                      {input.vcf}
        """)

rule _sv_merge_prepare_survivor:
    input:
        gridss_roi_filtered = str(rules._sv_merge_filter_gridss.output.vcf_roi_filtered),
        manta_roi_filtered = str(rules._sv_merge_filter_manta.output.vcf_roi_filtered)
    output:
        sample_sheet = CFG["dirs"]["survivor"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.survivor_sample_sheet.txt"
    shell:
        op.as_one_line("""
        ls {input.gridss_roi_filtered} {input.manta_roi_filtered} > {output.sample_sheet}
        """)
        
rule _sv_merge_run_survivor:
    input:
        sample_sheet = str(rules._sv_merge_prepare_survivor.output.sample_sheet),
    output:
        vcf = CFG["dirs"]["survivor"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.survivor.vcf"
    conda:
        CFG["conda_envs"]["survivor"]
    threads:
        CFG["threads"]["_sv_merge_run"]
    resources:
        **CFG["resources"]["_sv_merge_run"]
    shell:
        op.as_one_line("""
        SURVIVOR merge {input.sample_sheet} 1000 2 1 1 0 5 {output.vcf}
        """)
    
rule _sv_merge_format:
    input:
        gridss_roi_filtered = str(rules._sv_merge_filter_gridss.output.vcf_roi_filtered),
        manta_roi_filtered = str(rules._sv_merge_filter_manta.output.vcf_roi_filtered),
        survivor_vcf = str(rules._sv_merge_run_survivor.output.vcf),
        igcaller_tsv = str(rules._sv_merge_inputs.output.igcaller_tsv),
        roi_annotations = CFG["inputs"]["roi_annotations"],
    output:
        csv = CFG["dirs"]["merge_format"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_sv_merge_output.csv"
    params:
        script = CFG["scripts"]["format_output"]
    log:
        stdout = CFG["logs"]["merge_format"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.outputs.stdout.log",
        stderr = CFG["logs"]["merge_format"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.outputs.stderr.log"
    conda:
        CFG["conda_envs"]["postprocessing"]
    resources:
        **CFG["resources"]["_sv_merge_run"]
    shell:
        op.as_one_line("""
        Rscript --vanilla
                {params.script}
                -s {input.survivor_vcf}
                -g {input.gridss_roi_filtered}
                -m {input.manta_roi_filtered}
                -i {input.igcaller_tsv}
                -r {input.roi_annotations}
                -o {output.csv}
                > {log.stdout} 2> {log.stderr}
        """)
        
# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _sv_merge_output_tables:
    input:
        csv = str(rules._sv_merge_format.output.csv)
    output:
        csv = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_sv_merge_output.csv"
    run:
        op.relative_symlink(input.csv, output.csv, in_module= True)


rule _sv_merge_all:
    input:
        expand(
            [
                str(rules._sv_merge_output_tables.output.csv),
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
