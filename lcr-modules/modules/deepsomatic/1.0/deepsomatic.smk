#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Google/DeepSomatic Team
# Module Author:    Joanna A. Krupka
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["deepsomatic"]`
CFG = op.setup_module(
    name = "deepsomatic",
    version = "1.0",
    subdirectories = ["inputs", "deepsomatic", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _deepsomatic_input_bam,
    _deepsomatic_output_vcf,
    _deepsomatic_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _deepsomatic_input_bam:
    input:
        tumor_bam = CFG["inputs"]["sample_bam"],
        tumor_bai = CFG["inputs"]["sample_bai"],
        normal_bam = CFG["inputs"]["normal_bam"],
        normal_bai = CFG["inputs"]["normal_bai"]
    output:
        tumor_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.tumor.bam",
        tumor_bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.tumor.bam.bai",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.normal.bam",
        normal_bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.normal.bam.bai"
    run:
        op.relative_symlink(input.tumor_bam, output.tumor_bam)
        op.relative_symlink(input.tumor_bai, output.tumor_bai)
        # Handle optional normal BAM for tumor-only mode
        if input.normal_bam != "":
            op.relative_symlink(input.normal_bam, output.normal_bam)
            op.relative_symlink(input.normal_bai, output.normal_bai)
        else:
            # Create empty files for tumor-only mode
            shell("touch {output.normal_bam} {output.normal_bai}")


rule _deepsomatic_input_bam_tumor_only:
    input:
        tumor_bam = CFG["inputs"]["sample_bam"],
        tumor_bai = CFG["inputs"]["sample_bai"]
    output:
        tumor_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.tumor_only.bam",
        tumor_bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.tumor_only.bam.bai"
    run:
        op.relative_symlink(input.tumor_bam, output.tumor_bam)
        op.relative_symlink(input.tumor_bai, output.tumor_bai)


# Run DeepSomatic in matched tumor-normal mode
rule _deepsomatic_run_matched:
    input:
        tumor_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.tumor.bam",
        tumor_bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.tumor.bam.bai",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.normal.bam",
        normal_bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.normal.bam.bai",
        reference = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        reference_index = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai")
    output:
        vcf = CFG["dirs"]["deepsomatic"] + "{seq_type}--{genome_build}/{sample_id}.matched.deepsomatic.vcf.gz",
        tbi = CFG["dirs"]["deepsomatic"] + "{seq_type}--{genome_build}/{sample_id}.matched.deepsomatic.vcf.gz.tbi"
    params:
        output_dir = CFG["dirs"]["deepsomatic"] + "{seq_type}--{genome_build}/",
        sample_name_tumor = "{sample_id}_tumor",
        sample_name_normal = "{sample_id}_normal",
        model_type = lambda wildcards: CFG["options"]["model_type_matched"].get(wildcards.seq_type, "FFPE_WES"),
        regions = lambda wildcards: "--regions=" + CFG["options"]["regions"].get(wildcards.seq_type, "") if CFG["options"]["regions"].get(wildcards.seq_type, "") else "",
        extra_args = CFG["options"]["deepsomatic_matched"]
    container:
        CFG["options"]["docker_image"]
    threads:
        CFG["threads"]["deepsomatic_run"]
    resources:
        **CFG["resources"]["deepsomatic_run"]
    shell:
        """
        # Create output directory
        mkdir -p {params.output_dir}/intermediate_results_dir/{wildcards.sample_id}
        mkdir -p {params.output_dir}/logs/{wildcards.sample_id}

        # Run DeepSomatic in matched mode
        run_deepsomatic \
            --model_type={params.model_type} \
            --ref={input.reference} \
            --reads_normal={input.normal_bam} \
            --reads_tumor={input.tumor_bam} \
            --output_vcf={output.vcf} \
            --sample_name_tumor={params.sample_name_tumor} \
            --sample_name_normal={params.sample_name_normal} \
            --num_shards={threads} \
            --logging_dir={params.output_dir}/logs/{wildcards.sample_id} \
            --intermediate_results_dir={params.output_dir}/intermediate_results_dir/{wildcards.sample_id} \
            {params.regions} \
            {params.extra_args}
        """


# Run DeepSomatic in tumor-only mode
rule _deepsomatic_run_tumor_only:
    input:
        tumor_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.tumor_only.bam",
        tumor_bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.tumor_only.bam.bai",
        reference = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        reference_index = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai")
    output:
        vcf = CFG["dirs"]["deepsomatic"] + "{seq_type}--{genome_build}/{sample_id}.tumor_only.deepsomatic.vcf.gz",
        tbi = CFG["dirs"]["deepsomatic"] + "{seq_type}--{genome_build}/{sample_id}.tumor_only.deepsomatic.vcf.gz.tbi"
    params:
        output_dir = CFG["dirs"]["deepsomatic"] + "{seq_type}--{genome_build}/",
        sample_name_tumor = "{sample_id}_tumor",
        model_type = lambda wildcards: CFG["options"]["model_type_tumor_only"].get(wildcards.seq_type, "FFPE_WES_TUMOR_ONLY"),
        regions = lambda wildcards: "--regions=" + CFG["options"]["regions"].get(wildcards.seq_type, "") if CFG["options"]["regions"].get(wildcards.seq_type, "") else "",
        use_default_pon = CFG["options"]["use_default_pon_filtering"],
        extra_args = CFG["options"]["deepsomatic_tumor_only"]
    container:
        CFG["options"]["docker_image"]
    threads:
        CFG["threads"]["deepsomatic_run"]
    resources:
        **CFG["resources"]["deepsomatic_run"]
    shell:
        """
        # Create output directory
        mkdir -p {params.output_dir}/intermediate_results_dir/{wildcards.sample_id}
        mkdir -p {params.output_dir}/logs/{wildcards.sample_id}

        # Run DeepSomatic in tumor-only mode
        run_deepsomatic \
            --model_type={params.model_type} \
            --ref={input.reference} \
            --reads_tumor={input.tumor_bam} \
            --output_vcf={output.vcf} \
            --sample_name_tumor={params.sample_name_tumor} \
            --num_shards={threads} \
            --logging_dir={params.output_dir}/logs/{wildcards.sample_id} \
            --intermediate_results_dir={params.output_dir}/intermediate_results_dir/{wildcards.sample_id} \
            --use_default_pon_filtering={params.use_default_pon} \
            {params.regions} \
            {params.extra_args}
        """


# Symlink the final output files
rule _deepsomatic_output_vcf:
    input:
        vcf = CFG["dirs"]["deepsomatic"] + "{seq_type}--{genome_build}/{sample_id}.{mode}.deepsomatic.vcf.gz",
        tbi = CFG["dirs"]["deepsomatic"] + "{seq_type}--{genome_build}/{sample_id}.{mode}.deepsomatic.vcf.gz.tbi"
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{sample_id}.{mode}.deepsomatic.vcf.gz",
        tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{sample_id}.{mode}.deepsomatic.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf)
        op.relative_symlink(input.tbi, output.tbi)


# Generates the target sentinels for each module
rule _deepsomatic_all:
    input:
        expand(
            [
                CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{sample_id}.{mode}.deepsomatic.vcf.gz",
                CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{sample_id}.{mode}.deepsomatic.vcf.gz.tbi"
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"],
            mode=CFG["samples"]["mode"]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)