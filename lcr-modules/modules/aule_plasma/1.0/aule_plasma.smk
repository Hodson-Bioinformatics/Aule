#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Joanna A. Krupka, 
# Module Author:    Joanna A. Krupka, Ilias Moutsopoulos


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
# `CFG` is a shortcut to `config["lcr-modules"]["aule_plasma"]`
CFG = op.setup_module(
    name = "aule_plasma",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "sortbam","candidate_vcf", "annotated_vcf","outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _aule_plasma_input_vcf,
    _aule_plasma_input_bam_raw,
    _aule_plasma_input_bam_collapsed,
    _aule_plasma_candidate_variants,
    _aule_plasma_output_vcf,
    _aule_plasma_normalise, 
    _aule_plasma_annotate_pon,
    _aule_plasma_annotate_gnomad,
    _aule_plasma_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _aule_plasma_input_vcf:
    input:
        vcf = CFG["inputs"]["sample_vcf"],
        tbi = CFG["inputs"]["sample_tbi"],
    output:
        vcf = CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{base_name}.vcf.gz",
        tbi = CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{base_name}.vcf.tbi"
    run:
        op.absolute_symlink(input.vcf, output.vcf)
        op.absolute_symlink(input.tbi, output.tbi)

rule _aule_plasma_input_bam_raw:
    input:
        bam = CFG["inputs"]["sample_bam_raw"],
        bai = CFG["inputs"]["sample_bai_raw"],
    output:
        bam = CFG["dirs"]["inputs"] + "bam_raw/{seq_type}--{genome_build}/{sample_id}.raw.bam",
        bai = CFG["dirs"]["inputs"] + "bam_raw/{seq_type}--{genome_build}/{sample_id}.raw.bam.bai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        
rule _aule_plasma_input_bam_collapsed:
    input:
        bam = CFG["inputs"]["sample_bam_collapsed"],
    output:
        bam = CFG["dirs"]["inputs"] + "bam_collapsed/{seq_type}--{genome_build}/{sample_id}.collapsed.bam", 
    run:
        op.absolute_symlink(input.bam, output.bam)

#### Sort by coordinate 
rule _aule_plasma_sort_and_index:
    input:
        bam = str(rules._aule_plasma_input_bam_collapsed.output.bam)
    output:
        bam = temp(CFG["dirs"]["sortbam"] + "{seq_type}--{genome_build}/{sample_id}.sorted.bam"),
    log:
        stdout = CFG["logs"]["sortbam"] + "{seq_type}--{genome_build}/{sample_id}.sort_coordinates.log",
        stderr = CFG["logs"]["sortbam"] + "{seq_type}--{genome_build}/{sample_id}.sort_coordinates.stderr.log"
    params:
        tempdir = CFG["scratch_directory"]
    resources: 
        **CFG["resources"]["sortbam_coordinates"]
    threads:
        CFG["threads"]["sortbam_coordinates"]
    conda:
        CFG["conda_envs"]["picard"]
    shell:
        op.as_one_line("""
        picard SortSam {resources.mem_mb_java}
        --TMP_DIR {params.tempdir}
        -I {input.bam}
        -O {output.bam}
        -SO coordinate
        --CREATE_INDEX true
        > {log.stdout} 2> {log.stderr}
        """)
    
### Clean an input vcf file, remove all genotypes, keep candidate variants only 
rule _aule_plasma_candidate_variants:
    input:
        vcf = str(rules._aule_plasma_input_vcf.output.vcf)
    output:
        vcf = CFG["dirs"]["candidate_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{base_name}.variants_only.vcf.gz",
    log:
        log = CFG["logs"]["candidate_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{base_name}.variants_only.log"
    params:
        tempdir = CFG["scratch_directory"]
    resources: 
        **CFG["resources"]["bcftools"]
    threads:
        CFG["threads"]["bcftools"]
    conda:
        CFG["conda_envs"]["bcftools"]
    shell:
        op.as_one_line("""
        bcftools view -G -Oz {input.vcf} | bcftools annotate -x INFO,FORMAT -Oz | bcftools norm --atomize -m-both --threads {threads} -o {output.vcf} &&
        sleep 10 && 
        tabix -p vcf {output.vcf}
        """)
        
### Annotate the candidate variants with 
rule _aule_plasma_annotate_vcf:
    input:
        vcf = str(rules._aule_plasma_candidate_variants.output.vcf),
        tumour_bam = CFG["dirs"]["sortbam"] + "{seq_type}--{genome_build}/{tumour_id}.sorted.bam",
        normal_bam = CFG["dirs"]["sortbam"] + "{seq_type}--{genome_build}/{normal_id}.sorted.bam"
    output:
        vcf = CFG["dirs"]["annotated_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{base_name}.scored.vcf.gz",
        tbi = CFG["dirs"]["annotated_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{base_name}.scored.vcf.gz.tbi",
    log:
        log = CFG["logs"]["annotated_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{base_name}.scored.log"
    params:
        script = CFG["scripts"]["annotate_variants"],
        genome_build = "{genome_build}",
        min_base_qual = CFG["options"]["min_base_qual"],
        min_mapping_qual = CFG["options"]["min_mapping_qual"]
    resources: 
        **CFG["resources"]["annotated_variants"]
    threads:
        CFG["threads"]["annotated_variants"]
    conda:
        CFG["conda_envs"]["pysam"]
    shell:
        op.as_one_line("""
        python {params.script}
        -vcf_in {input.vcf} 
        -vcf_out {output.vcf} 
        -tumour_bam {input.tumour_bam} 
        -normal_bam {input.normal_bam} 
        -threads {threads}
        -genome_build {params.genome_build} 
        -min_base_qual {params.min_base_qual} 
        -min_mapping_qual {params.min_mapping_qual} && 
        tabix -p vcf {output.vcf}
        """)
        
rule _aule_plasma_annotate_gnomad:
    input:
        vcf = str(rules._aule_plasma_annotate_vcf.output.vcf),
        normalized_gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.normalized.{genome_build}.vcf.gz")
    output:
        vcf = CFG["dirs"]["annotated_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{base_name}.gnomad.vcf.gz",
        tbi = CFG["dirs"]["annotated_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{base_name}.gnomad.vcf.gz.tbi",
    conda:
        CFG["conda_envs"]["bcftools"]
    resources:
        **CFG["resources"]["bcftools"]
    threads:
        CFG["threads"]["bcftools"]
    shell:
        op.as_one_line("""
        bcftools annotate --threads {threads} -a {input.normalized_gnomad} {input.vcf} -c "INFO/gnomADg_AF:=INFO/AF" -o {output.vcf} && 
        tabix -p vcf {output.vcf}
        """)
        
rule _aule_plasma_annotate_pon:
    input:
        vcf = str(rules._aule_plasma_annotate_gnomad.output.vcf),
        pon = CFG["inputs"]["pon"]
    output:
        vcf = CFG["dirs"]["annotated_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{base_name}.pon.vcf.gz",
        tbi = CFG["dirs"]["annotated_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{base_name}.pon.vcf.gz.tbi",
    conda:
        CFG["conda_envs"]["bcftools"]
    resources:
        **CFG["resources"]["bcftools"]
    threads:
        CFG["threads"]["bcftools"]
    shell:
        op.as_one_line("""
        bcftools annotate --threads {threads} -a {input.pon} {input.vcf} -c "INFO/pon_AF:=INFO/FRACTION" -o {output.vcf} && 
        tabix -p vcf {output.vcf}
        """)
   
# rule _aule_plasma_filter: 
#     input:
#         vcf = str(rules._aule_plasma_candidate_variants.output.vcf),
#     output:
#         vcf = CFG["dirs"]["annotated_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{base_name}.scored.vcf.gz",
#         tbi = CFG["dirs"]["annotated_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{base_name}.scored.vcf.gz.tbi",
#     log:
#         log = CFG["logs"]["annotated_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{base_name}.scored.log"
#     params:
#         script = CFG["scripts"]["annotate_variants"],
#         genome_build = "{genome_build}",
#         min_base_qual = CFG["options"]["min_base_qual"],
#         min_mapping_qual = CFG["options"]["min_mapping_qual"]
#     resources: 
#         **CFG["resources"]["annotated_variants"]
#     threads:
#         CFG["threads"]["annotated_variants"]
#     conda:
#         CFG["conda_envs"]["pysam"]
#     shell:
#         op.as_one_line("""
#         python {params.script}
#         -vcf_in {input.vcf} 
#         -vcf_out {output.vcf} 
#         -tumour_bam {input.tumour_bam} 
#         -normal_bam {input.normal_bam} 
#         -threads {threads}
#         -genome_build {params.genome_build} 
#         -min_base_qual {params.min_base_qual} 
#         -min_mapping_qual {params.min_mapping_qual} && 
#         tabix -p vcf {output.vcf}
#         """)

# # Symlinks the final output files into the module results directory (under '99-outputs/')
rule _aule_plasma_output_vcf:
    input:
        vcf = str(rules._aule_plasma_annotate_pon.output.vcf),
        tbi = str(rules._aule_plasma_annotate_pon.output.tbi)
    output:
        vcf = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{base_name}.low_vaf.vcf.gz",
        tbi = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{base_name}.low_vaf.vcf.gz.tbi",
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module= True)
        op.relative_symlink(input.tbi, output.tbi, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _aule_plasma_all:
    input:
        expand(
            [
                str(rules._aule_plasma_output_vcf.output.vcf),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"],
            base_name = [CFG["vcf_base_name"]] * len(CFG["runs"]["tumour_sample_id"]),
            allow_missing = True)


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
