#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Joanna A. Krupka
# Module Author:    Joanna A. Krupka
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
# `CFG` is a shortcut to `config["lcr-modules"]["aule_snv_pon"]`
CFG = op.setup_module(
    name = "aule_snv_pon",
    version = "1.0",
    subdirectories = ["inputs", "mutect2", "genomicsdb", "pon", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _aule_snv_pon_input_bam,
    _aule_snv_pon_output_vcf,
    _aule_snv_pon_all


##### RULES #####
def _mutect_get_capspace(wildcards):
    CFG = config["lcr-modules"]["aule_snv_pon"]
    try:
        # Get the appropriate capture space for this sample
        cap_space = op.get_capture_space(CFG, wildcards.sample_id, wildcards.genome_build, wildcards.seq_type, "interval_list")
        cap_space = reference_files(cap_space)
        this_space = cap_space
    # If we are using an older version of the reference workflow, we don't need to do anything
    except NameError:
        this_space = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.bed")
    # If this isn't a capture sample, we don't have a capture space, so return nothing
    if wildcards.seq_type != "capture":
        this_space = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.bed")
    return this_space

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _aule_snv_pon_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)

# Run Mutect2 in tumour-only mode
rule _aule_snv_pon_mutect2:
    input:
        normal_bam = str(rules._aule_snv_pon_input_bam.output.bam),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        dict = reference_files("genomes/{genome_build}/genome_fasta/genome.dict"),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz"),
        capture_arg = _mutect_get_capspace
    output:
        vcf = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}/output.vcf.gz",
        tbi = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}/output.vcf.gz.tbi",
        stat = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}/output.vcf.gz.stats"
    log:
        stdout = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}/mutect2_run.stdout.log",
        stderr = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}/mutect2_run.stderr.log"
    resources:
        **CFG["resources"]["mutect2_run"]
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        opts = CFG["options"]["mutect2_run"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["mutect2_run"]
    wildcard_constraints: 
        pair_status = "no_normal"
    shell:
        op.as_one_line("""
        gatk Mutect2 --java-options "-Xmx{params.mem_mb}m" 
        {params.opts} -I {input.normal_bam} -R {input.fasta} 
        -O {output.vcf} --germline-resource {input.gnomad} 
        -L {input.capture_arg}
        > {log.stdout} 2> {log.stderr}
        """)

# Generate PON genomicsDB
def _aule_snv_pon_get_vcfs(wildcards):
    CFG = config["lcr-modules"]["aule_snv_pon"]
    samples = CFG["samples"]["sample_id"]
    vcfs = expand(
        [CFG["dirs"]["mutect2"] + "{{seq_type}}--{{genome_build}}/{sample_id}/output.vcf.gz"],
        sample_id = samples
    )
    return(vcfs)

        
rule _aule_snv_pon_genomicsdb: 
    input: 
        vcfs = _aule_snv_pon_get_vcfs,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
    output: 
        genomicsdb = directory(CFG["dirs"]["genomicsdb"] + "{seq_type}--{genome_build}/pon_db")
    log:
        stdout = CFG["logs"]["genomicsdb"] + "{seq_type}--{genome_build}/genomicsdb.stdout.log",
        stderr = CFG["logs"]["genomicsdb"] + "{seq_type}--{genome_build}/genomicsdb.stderr.log"
    resources:
        **CFG["resources"]["genomicsdb"]
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        opts = CFG["options"]["genomicsdb"],
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["genomicsdb"]
    shell:
        op.as_one_line("""
        gatk GenomicsDBImport --java-options "-Xmx{params.mem_mb}m" 
        {params.opts} 
        $(for input in {input.vcfs}; do printf -- "-V $input "; done)
        --genomicsdb-workspace-path {output.genomicsdb}
        > {log.stdout} 2> {log.stderr}
        """)
        
rule _aule_snv_pon_createpon:        
    input: 
        genomicsdb = str(rules._aule_snv_pon_genomicsdb.output.genomicsdb),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output: 
        pon = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/pon.vcf.gz"
    log:
        stdout = CFG["logs"]["pon"] + "{seq_type}--{genome_build}/createpon.stdout.log",
        stderr = CFG["logs"]["pon"] + "{seq_type}--{genome_build}/createpon.stderr.log"
    resources:
        **CFG["resources"]["createpon"]
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        opts = CFG["options"]["createpon"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["createpon"]
    shell:
        op.as_one_line("""
        gatk CreateSomaticPanelOfNormals --java-options "-Xmx{params.mem_mb}m" 
        {params.opts} 
        -V gendb://{input.genomicsdb}
        -R {input.fasta}
        -O {output.pon}
        > {log.stdout} 2> {log.stderr}
        """)
        

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _aule_snv_pon_output_vcf:
    input:
        pon = str(rules._aule_snv_pon_createpon.output.pon)
    output:
        pon = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/pon.vcf.gz"
    run:
        op.relative_symlink(input.pon, output.pon, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _aule_snv_pon_all:
    input:
        expand(
            [
                str(rules._aule_snv_pon_output_vcf.output.pon)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
