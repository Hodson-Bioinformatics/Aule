#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ilias Moutsopoulos 
# Module Author:    Ilias Moutsopoulos
# Contributors:     Joanna Krupka


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
    version = "1.1",
    subdirectories = ["inputs", "sort_bam", "annotate_vcf_pon_error", "annotate_vcf_slms3_umi", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _aule_plasma_input_vcf,
    _aule_plasma_input_bam,
    _aule_plasma_input_positional_error_rates,
    _aule_plasma_input_directional_error_rates,
    _aule_plasma_output_vcf,
    _aule_plasma_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _aule_plasma_input_vcf:
    input:
        vcf_2plus = CFG["inputs"]["sample_vcf_2plus"],
        tbi_2plus = CFG["inputs"]["sample_tbi_2plus"],
        vcf_3plus = CFG["inputs"]["sample_vcf_3plus"],
        tbi_3plus = CFG["inputs"]["sample_tbi_3plus"]
    output:
        vcf_2plus = CFG["dirs"]["inputs"] + "vcf/2plus/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.2+.vcf.gz",
        tbi_2plus = CFG["dirs"]["inputs"] + "vcf/2plus/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.2+.vcf.gz.tbi",
        vcf_3plus = CFG["dirs"]["inputs"] + "vcf/3plus/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.3+.vcf.gz",
        tbi_3plus = CFG["dirs"]["inputs"] + "vcf/3plus/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.3+.vcf.gz.tbi"
    run:
        op.absolute_symlink(input.vcf_2plus, output.vcf_2plus)
        op.absolute_symlink(input.tbi_2plus, output.tbi_2plus)
        op.absolute_symlink(input.vcf_3plus, output.vcf_3plus)
        op.absolute_symlink(input.tbi_3plus, output.tbi_3plus)

rule _aule_plasma_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam"
    run:
        op.absolute_symlink(input.bam, output.bam)

rule _aule_plasma_input_pon_vcf:
    input:
        pon = CFG["inputs"]["pon_vcf"],
        tbi = CFG["inputs"]["pon_tbi"]
    output:
        pon = CFG["dirs"]["inputs"] + "pon/{seq_type}--{genome_build}/pon.vcf.gz",
        tbi = CFG["dirs"]["inputs"] + "pon/{seq_type}--{genome_build}/pon.vcf.gz.tbi",
    run:
        op.absolute_symlink(input.pon, output.pon)
        op.absolute_symlink(input.tbi, output.tbi)
        
rule _aule_plasma_input_positional_error_rates:
    input:
        positional_error_rates = CFG["inputs"]["positional_error_rates"],
        positional_error_header = CFG["inputs"]["positional_error_header"],
        positional_error_tbi = CFG["inputs"]["positional_error_tbi"]
    output:
        positional_error_rates = CFG["dirs"]["inputs"] + "error_rates/{seq_type}--{genome_build}/positional_error_rates.txt.gz",
        positional_error_header = CFG["dirs"]["inputs"] + "error_rates/{seq_type}--{genome_build}/positional_error_header.txt",
        positional_error_tbi = CFG["dirs"]["inputs"] + "error_rates/{seq_type}--{genome_build}/positional_error_rates.txt.gz.tbi",
    run:
        op.absolute_symlink(input.positional_error_rates, output.positional_error_rates)
        op.absolute_symlink(input.positional_error_header, output.positional_error_header)
        op.absolute_symlink(input.positional_error_tbi, output.positional_error_tbi)

rule _aule_plasma_input_directional_error_rates:
    input:
        directional_error_rates = CFG["inputs"]["directional_error_rates"],
        directional_error_header = CFG["inputs"]["directional_error_header"],
        directional_error_tbi = CFG["inputs"]["directional_error_tbi"],
    output:
        directional_error_rates = CFG["dirs"]["inputs"] + "error_rates/{seq_type}--{genome_build}/directional_error_rates.txt.gz",
        directional_error_header = CFG["dirs"]["inputs"] + "error_rates/{seq_type}--{genome_build}/directional_error_header.txt",
        directional_error_tbi = CFG["dirs"]["inputs"] + "error_rates/{seq_type}--{genome_build}/directional_error_rates.txt.gz.tbi",
    run:
        op.absolute_symlink(input.directional_error_rates, output.directional_error_rates)
        op.absolute_symlink(input.directional_error_header, output.directional_error_header)
        op.absolute_symlink(input.directional_error_tbi, output.directional_error_tbi)


# Sort bam by coordinate 
rule _aule_plasma_sort_bam:
    input:
        bam = str(rules._aule_plasma_input_bam.output.bam)
    output:
        bam = temp(CFG["dirs"]["sort_bam"] + "{seq_type}--{genome_build}/{tumour_id}.sorted.bam"),
    log:
        stdout = CFG["logs"]["sort_bam"] + "{seq_type}--{genome_build}/{tumour_id}.sort_bam.stdout.log",
        stderr = CFG["logs"]["sort_bam"] + "{seq_type}--{genome_build}/{tumour_id}.sort_bam.stderr.log"
    params:
        tempdir = CFG["scratch_directory"]
    resources: 
        **CFG["resources"]["sort_bam"]
    threads:
        CFG["threads"]["sort_bam"]
    conda:
        CFG["conda_envs"]["picard"]
    shell:
        op.as_one_line("""
       if samtools view -H {input.bam} | grep -q 'SO:coordinate'; then 
            ln -sf ../../../../{input.bam} {output.bam} && sleep 30 && samtools index {output.bam};
        else 
            picard SortSam 
                --TMP_DIR {params.tempdir} 
                -I {input.bam} 
                -O {output.bam} 
                -SO coordinate 
                --CREATE_INDEX true ;
        fi 
        > {log.stdout} 2> {log.stderr}
        """)

# Annotate VCF with PON fractions 
rule _aule_plasma_annotate_vcf_pon_error:
    input:
        vcf_2plus = str(rules._aule_plasma_input_vcf.output.vcf_2plus),
        directional_error_rates = str(rules._aule_plasma_input_directional_error_rates.output.directional_error_rates),
        directional_error_header = str(rules._aule_plasma_input_directional_error_rates.output.directional_error_header),
        positional_error_rates = str(rules._aule_plasma_input_positional_error_rates.output.positional_error_rates),
        positional_error_header = str(rules._aule_plasma_input_positional_error_rates.output.positional_error_header),
        pon = str(rules._aule_plasma_input_pon_vcf.output.pon),
    output:
        vcf_2plus = CFG["dirs"]["annotate_vcf_pon_error"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.annotate_pon_error.vcf.gz",
        tbi = CFG["dirs"]["annotate_vcf_pon_error"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.annotate_pon_error.vcf.gz.tbi",
    log:
        log = CFG["logs"]["annotate_vcf_pon_error"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.annotate_vcf_pon_error.log"
    params:
        pon_af_field = CFG["options"]["pon_af_field"]
    resources: 
        **CFG["resources"]["annotate_vcf_pon_error"]
    threads:
        CFG["threads"]["annotate_vcf_pon_error"]
    conda:
        CFG["conda_envs"]["bcftools"]
    shell:
        op.as_one_line("""
        bcftools annotate -a {input.pon} -c {params.pon_af_field} {input.vcf_2plus} |
        bcftools annotate -c INFO/PON_AF:=INFO/{params.pon_af_field} - |
        bcftools annotate  -a {input.positional_error_rates} -h {input.positional_error_header} -c CHROM,POS,umi_error_rate,r12_rate - |
        bcftools annotate  -a {input.directional_error_rates} -h {input.directional_error_header} -c CHROM,POS,REF,ALT,directional_error_rate - -Oz -o {output.vcf_2plus} &&
        sleep 30 &&
        tabix -p vcf {output.vcf_2plus}
        """)

# Annotate the variants in the 2+ vcf with UMI information, as well as whether they are in the 3+ vcf 
rule _aule_plasma_annotate_vcf_slms3_umi:
    input:
        vcf_2plus = str(rules._aule_plasma_annotate_vcf_pon_error.output.vcf_2plus),
        vcf_3plus = str(rules._aule_plasma_input_vcf.output.vcf_3plus),
        bam = str(rules._aule_plasma_sort_bam.output.bam)
    output:
        vcf = CFG["dirs"]["annotate_vcf_slms3_umi"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.annotate_slms3_umi.vcf.gz",
        tbi = CFG["dirs"]["annotate_vcf_slms3_umi"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.annotate_slms3_umi.vcf.gz.tbi"
    log:
        log = CFG["logs"]["annotate_vcf_slms3_umi"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.annotate_vcf_slms3_umi.log"
    params:
        script = CFG["scripts"]["annotate_vcf_slms3_umi"],
    resources: 
        **CFG["resources"]["annotate_vcf_slms3_umi"]
    threads:
        CFG["threads"]["annotate_vcf_slms3_umi"]
    conda:
        CFG["conda_envs"]["pysam"]
    shell:
        op.as_one_line("""
        python {params.script}
            -input_vcf_path {input.vcf_2plus} 
            -bam_file {input.bam} 
            -vcf_3plus {input.vcf_3plus} 
            -output_vcf_path {output.vcf} && 
        tabix -p vcf {output.vcf}
        """)


# Create a symlink to the final output VCF file
rule _aule_plasma_output_vcf:
    input:
        vcf = str(rules._aule_plasma_annotate_vcf_slms3_umi.output.vcf),
        tbi = str(rules._aule_plasma_annotate_vcf_slms3_umi.output.tbi)
    output:
        vcf = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.annotated.vcf.gz",
        tbi = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.annotated.vcf.gz.tbi"
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
            pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
