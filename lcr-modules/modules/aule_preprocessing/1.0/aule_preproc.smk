#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Joanna A. Krupka
# Module Author:    Joanna A. Krupka
# Contributors:     Ilias Moutsopoulos


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
# `CFG` is a shortcut to `config["lcr-modules"]["aule_preprocessing"]`
CFG = op.setup_module(
    name = "aule_preprocessing",
    version = "1.0",
    subdirectories = ["inputs", "fastq_to_unmapped_bam", "annotate_bam_with_umis", "mark_illumina_adapters", "bwa_input", "bwa_mem", "merge_bam_alignment", "fix_tags", "mark_duplicates", "base_recalibration", "apply_bsqr", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _aule_preprocessing_input_fastq,
    _aule_preprocessing_output_bam,
    _aule_preprocessing_all

##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _aule_preprocessing_input_fastq:
    input:
        fastq_R1 = lambda wildcards: SAMPLES.loc[SAMPLES["sample_id"] == wildcards.sample_id, "sample_fastq_R1"],
        fastq_R2 = lambda wildcards: SAMPLES.loc[SAMPLES["sample_id"] == wildcards.sample_id, "sample_fastq_R2"],
        fastq_UMI = lambda wildcards: SAMPLES.loc[SAMPLES["sample_id"] == wildcards.sample_id, "sample_fastq_UMI"]
    output:
        fastq_R1 = CFG["dirs"]["inputs"] + "fastq/{seq_type}/{sample_id}_R1.fastq.gz",
        fastq_R2 = CFG["dirs"]["inputs"] + "fastq/{seq_type}/{sample_id}_R2.fastq.gz",
        fastq_UMI = CFG["dirs"]["inputs"] + "fastq/{seq_type}/{sample_id}_UMI.fastq.gz"
    run:
        op.absolute_symlink(input.fastq_R1, output.fastq_R1)
        op.absolute_symlink(input.fastq_R2, output.fastq_R2)
        op.absolute_symlink(input.fastq_UMI, output.fastq_UMI)

# Converting FASTQ files to uBAMs with sample metadata
rule _aule_preprocessing_fastq_to_unmapped_bam:
    priority: -8
    input:
        fastq_R1 = str(rules._aule_preprocessing_input_fastq.output.fastq_R1),
        fastq_R2 = str(rules._aule_preprocessing_input_fastq.output.fastq_R2)
    output:
        ubam = temp(CFG["dirs"]["fastq_to_unmapped_bam"] + "{seq_type}/{sample_id}.unmapped.bam")
    log:
        stdout = CFG["logs"]["fastq_to_unmapped_bam"] + "{seq_type}/{sample_id}.preproc_fastq_to_unmapped.stdout.log",
        stderr = CFG["logs"]["fastq_to_unmapped_bam"] + "{seq_type}/{sample_id}.preproc_fastq_to_unmapped.stderr.log"
    params:
        options = CFG["options"]["fastq_to_unmapped_bam"],
        library_name = lambda wildcards: SAMPLES.loc[SAMPLES["sample_id"] == wildcards.sample_id, "CMDL_ID"].values[0],
        comment = lambda wildcards: SAMPLES.loc[SAMPLES["sample_id"] == wildcards.sample_id, "seq_batch"].values[0],
        scratch = CFG["scratch_directory"]     
    conda:
        CFG["conda_envs"]["picard"]
    threads:
        CFG["threads"]["fastq_to_unmapped_bam"]
    resources:
        **CFG["resources"]["fastq_to_unmapped_bam"]
    shell: 
        op.as_one_line(""" 
        picard FastqToSam 
            --FASTQ {input.fastq_R1} 
            --FASTQ2 {input.fastq_R2} 
            --OUTPUT {output.ubam} 
            --SAMPLE_NAME {wildcards.sample_id}
            --LIBRARY_NAME {params.library_name} 
            --COMMENT {params.comment} 
            --TMP_DIR {params.scratch}
            {params.options}
            > {log.stdout} 2> {log.stderr}
            """)
    
# Annotate UMIs to RX tag 
rule _aule_preprocessing_annotate_bam_with_umis: 
    priority: -7
    input:
        ubam = str(rules._aule_preprocessing_fastq_to_unmapped_bam.output.ubam),
        UMI = str(rules._aule_preprocessing_input_fastq.output.fastq_UMI)
    output:
        ubam = temp(CFG["dirs"]["annotate_bam_with_umis"] + "{seq_type}/{sample_id}.unmapped_umi.bam")
    log:
        stdout = CFG["logs"]["annotate_bam_with_umis"] + "{seq_type}/{sample_id}.annotate_bam_with_umis.stdout.log",
        stderr = CFG["logs"]["annotate_bam_with_umis"] + "{seq_type}/{sample_id}.annotate_bam_with_umis.stderr.log"
    params:
        options = CFG["options"]["annotate_bam_with_umis"],
        scratch = CFG["scratch_directory"],
        java = CFG["options"]["java_fgbio"]
    conda:
        CFG["conda_envs"]["fgbio"]
    threads:
        CFG["threads"]["annotate_bam_with_umis"]
    resources:
        **CFG["resources"]["annotate_bam_with_umis"]
    shell:
        op.as_one_line("""
        fgbio --tmp-dir {params.scratch} {params.java} 
        AnnotateBamWithUmis 
        -i {input.ubam} 
        -f {input.UMI} 
        -o {output.ubam} 
        {params.options} 
        > {log.stdout} 2> {log.stderr}
        """)

# Find adapters
rule _aule_preprocessing_mark_illumina_adapters:
    priority: -6
    input:
        ubam = str(rules._aule_preprocessing_annotate_bam_with_umis.output.ubam)
    output:
        ubam = temp(CFG["dirs"]["mark_illumina_adapters"] + "{seq_type}/{sample_id}.unmapped_umi_marked.bam")
    log:
        raport = CFG["logs"]["mark_illumina_adapters"] + "{seq_type}/{sample_id}.mark_illumina_adapters.txt",
        stdout = CFG["logs"]["mark_illumina_adapters"] + "{seq_type}/{sample_id}.mark_illumina_adapters.stdout.log",
        stderr = CFG["logs"]["mark_illumina_adapters"] + "{seq_type}/{sample_id}.mark_illumina_adapters.stderr.log"
    conda:
        CFG["conda_envs"]["picard"]
    threads:
        CFG["threads"]["mark_illumina_adapters"]
    resources:
        **CFG["resources"]["mark_illumina_adapters"]
    shell:
        op.as_one_line("""
        picard MarkIlluminaAdapters \
        --INPUT {input.ubam} \
        --METRICS {log.raport} \
        --OUTPUT {output.ubam}
         > {log.stdout} 2> {log.stderr}
        """)

# BWA input
rule _aule_preprocessing_bwa_input:
    priority: -5
    input:
        ubam = str(rules._aule_preprocessing_mark_illumina_adapters.output.ubam)
    output:
        fastq = temp(CFG["dirs"]["bwa_input"] + "{seq_type}/{sample_id}.bwa_input.fastq")
    log:
        stdout = CFG["logs"]["bwa_input"] + "{seq_type}/{sample_id}.bwa_input.stdout.log",
        stderr = CFG["logs"]["bwa_input"] + "{seq_type}/{sample_id}.bwa_input.stderr.log"
    params:
        options = CFG["options"]["bwa_input"]
    conda:
        CFG["conda_envs"]["picard"]
    threads:
        CFG["threads"]["bwa_input"]
    resources:
        **CFG["resources"]["bwa_input"]
    shell:
        op.as_one_line("""
        picard SamToFastq
        INPUT={input.ubam}
        FASTQ={output.fastq}
        {params.options}
         > {log.stdout} 2> {log.stderr}
        """)

# BWA MEM (piped to samtools view to save the output as BAM)
rule _aule_preprocessing_bwa_mem:
    priority: -4
    input:
        fastq = str(rules._aule_preprocessing_bwa_input.output.fastq),
        ref = reference_files("genomes/{genome_build}/bwa_index/bwa-0.7.17/genome.fa")
    output:
        bam = temp(CFG["dirs"]["bwa_mem"] + "{seq_type}--{genome_build}/{sample_id}.bwa.bam")
    log:
        stderr = CFG["logs"]["bwa_mem"] + "{seq_type}--{genome_build}/{sample_id}.bwa.stderr.log"
    params:
        options = CFG["options"]["bwa_mem"]
    conda:
        CFG["conda_envs"]["bwa_aule"] 
    threads:
        CFG["threads"]["bwa_mem"]
    resources:
        **CFG["resources"]["bwa_mem"]
    shell:
        op.as_one_line("""
        bwa mem {params.options}
        -t {threads} 
        {input.ref} 
        {input.fastq} | 
        samtools view -1 - > {output.bam}
        2> {log.stderr}
        """)

# Merge ubam and bam
rule _aule_preprocessing_merge_bam_alignment:
    priority: -3
    input:
        ubam = str(rules._aule_preprocessing_mark_illumina_adapters.output.ubam),
        bam = str(rules._aule_preprocessing_bwa_mem.output.bam),
        ref = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        bam = temp(CFG["dirs"]["merge_bam_alignment"] + "{seq_type}--{genome_build}/{sample_id}.merged.bam")
    log:
        stdout = CFG["logs"]["merge_bam_alignment"] + "{seq_type}--{genome_build}/{sample_id}.merge_bam_alignment.stdout.log",
        stderr = CFG["logs"]["merge_bam_alignment"] + "{seq_type}--{genome_build}/{sample_id}.merge_bam_alignment.stderr.log"
    params:
        options = CFG["options"]["merge_bam_alignment"],
        scratch = CFG["scratch_directory"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["merge_bam_alignment"]
    resources:
        **CFG["resources"]["merge_bam_alignment"]
    shell:
         op.as_one_line("""
        gatk MergeBamAlignment 
       --ALIGNED_BAM {input.bam} 
       --UNMAPPED_BAM {input.ubam} 
       --OUTPUT {output.bam} 
       --REFERENCE_SEQUENCE {input.ref} 
       {params.options} 
       --TMP_DIR {params.scratch}
       > {log.stdout} 2> {log.stderr}
       """)

# Fix tag values for NM and UQ 
rule _aule_preprocessing_fix_tags:
    priority: -2
    input:
        bam = str(rules._aule_preprocessing_merge_bam_alignment.output.bam),
        ref = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output: 
        bam = temp(CFG["dirs"]["fix_tags"] + "{seq_type}--{genome_build}/{sample_id}.fixed.bam")
    log:
        stdout = CFG["logs"]["fix_tags"] + "{seq_type}--{genome_build}/{sample_id}.fix_tags.stdout.log",
        stderr = CFG["logs"]["fix_tags"] + "{seq_type}--{genome_build}/{sample_id}.fix_tags.stderr.log"
    params:
        options = CFG["options"]["fix_tags"]
    conda:
        CFG["conda_envs"]["gatk"] 
    threads:
        CFG["threads"]["fix_tags"]
    resources:
        **CFG["resources"]["fix_tags"]
    shell:
        op.as_one_line("""
        gatk SetNmMdAndUqTags 
       --INPUT {input.bam}
       --OUTPUT {output.bam} 
       --REFERENCE_SEQUENCE {input.ref}
       {params.options}
        > {log.stdout} 2> {log.stderr}
        """)
        
# Mark Duplicates
rule _aule_preprocessing_mark_duplicates:
    priority: -1
    input:
        bam = str(rules._aule_preprocessing_fix_tags.output.bam),
    output:
        bam = temp(CFG["dirs"]["mark_duplicates"] + "{seq_type}--{genome_build}/{sample_id}.deduplicated.bam")
    log:
        metrics = CFG["logs"]["mark_duplicates"] + "{seq_type}--{genome_build}/{sample_id}.mark_duplicates.txt",
        stdout = CFG["logs"]["mark_duplicates"] + "{seq_type}--{genome_build}/{sample_id}.mark_duplicates.stdout.log",
        stderr = CFG["logs"]["mark_duplicates"] + "{seq_type}--{genome_build}/{sample_id}.mark_duplicates.stderr.log"
    params:
        options = CFG["options"]["mark_duplicates"],
        java = CFG["options"]["java_gatk"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["mark_duplicates"]
    resources:
        **CFG["resources"]["mark_duplicates"]
    shell:
        op.as_one_line("""
        gatk  {params.java} MarkDuplicates
       --INPUT {input.bam}
       --OUTPUT {output.bam}
       --METRICS_FILE {log.metrics}
       {params.options}
        > {log.stdout} 2> {log.stderr}
        """)

# Base recalibration
rule _aule_preprocessing_base_recalibraton:
    input:
        bam = str(rules._aule_preprocessing_mark_duplicates.output.bam),
        ref = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        region = lambda wildcards: config["lcr-modules"]["aule_preprocessing"]["inputs"]["bed"][wildcards.genome_build]
    output:
        recalibration = CFG["dirs"]["base_recalibration"] + "{seq_type}--{genome_build}/{sample_id}.base_recalibration.txt"
    log:
        stdout = CFG["logs"]["base_recalibration"] + "{seq_type}--{genome_build}/{sample_id}.base_recalibration.stdout.log",
        stderr = CFG["logs"]["base_recalibration"] + "{seq_type}--{genome_build}/{sample_id}.base_recalibration.stderr.log"
    params:
        options = CFG["options"]["base_recalibration"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["base_recalibration"]
    resources:
        **CFG["resources"]["base_recalibration"]
    shell:
        op.as_one_line("""
        gatk BaseRecalibrator 
        -R {input.ref} 
        -I {input.bam} 
        -L {input.region} 
        -O {output.recalibration}
        {params.options}
         > {log.stdout} 2> {log.stderr}
        """)

# Apply base recalibration
rule _aule_preprocessing_apply_bsqr:
    input:
        bam = str(rules._aule_preprocessing_mark_duplicates.output.bam),
        recalibration = str(rules._aule_preprocessing_base_recalibraton.output.recalibration),
        ref = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        bam = CFG["dirs"]["apply_bsqr"] + "{seq_type}--{genome_build}/{sample_id}.recalibrated.bam",
        bai = CFG["dirs"]["apply_bsqr"] + "{seq_type}--{genome_build}/{sample_id}.recalibrated.bai",
        md5 = CFG["dirs"]["apply_bsqr"] + "{seq_type}--{genome_build}/{sample_id}.recalibrated.bam.md5"
    log:
        stdout = CFG["logs"]["apply_bsqr"] + "{seq_type}--{genome_build}/{sample_id}.applyBSQR.stdout.log",
        stderr = CFG["logs"]["apply_bsqr"] + "{seq_type}--{genome_build}/{sample_id}.applyBSQR.stderr.log"
    params:
        options = CFG["options"]["apply_bsqr"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["apply_bsqr"]
    resources:
        **CFG["resources"]["apply_bsqr"]
    shell:
        op.as_one_line("""
        gatk ApplyBQSR 
        -R {input.ref} 
        -I {input.bam} 
        -O {output.bam} 
        -bqsr {input.recalibration} 
        {params.options}
         > {log.stdout} 2> {log.stderr}
        """)
        
# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _aule_preprocessing_output_bam:
    input:
        bam = str(rules._aule_preprocessing_apply_bsqr.output.bam),
        bai = str(rules._aule_preprocessing_apply_bsqr.output.bai),
        md5 = str(rules._aule_preprocessing_apply_bsqr.output.md5)
    output:
        bam = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.bam.bai",
        md5 = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.bam.md5",
    run:
        op.relative_symlink(input.bam, output.bam, in_module = True)
        op.relative_symlink(input.bai, output.bai, in_module = True)
        op.relative_symlink(input.md5, output.md5, in_module = True)

# Generates the target sentinels for each run, which generate the symlinks
rule _aule_preprocessing_all:
    input:
        expand(
            [
                str(rules._aule_preprocessing_output_bam.output.bam),
                str(rules._aule_preprocessing_output_bam.output.bai),
                str(rules._aule_preprocessing_output_bam.output.md5)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])
            
##### CLEANUP #####

# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
