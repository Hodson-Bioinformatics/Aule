#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  NA
# Module Author:    Prasath Pararajalingam
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
# `CFG` is a shortcut to `config["lcr-modules"]["igcaller"]`
CFG = op.setup_module(
    name = "igcaller",
    version = "1.0",
    subdirectories = ["inputs", "igcaller", "outputs"]
)


# Define rules to be run locally when using a compute cluster
localrules:
    _igcaller_input_bam,
    _igcaller_output_tsv,
    _igcaller_install,
    _igcaller_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _igcaller_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.crai"
    group: 
        "input_and_igcaller_run"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)


# Download and install IgCaller v1.3 if it doesn't exist in CFG["igcaller_path"]
rule _igcaller_install:
    params:
        igcaller = CFG["igcaller_path"]
    output:
        complete = CFG["igcaller_path"] + "igcaller-1.3.installed" # MODIFIED to v1.3
    shell:
        """
        if [ ! -f {params.igcaller}/IgCaller-1.3/IgCaller.py ]; then # MODIFIED to v1.3
            mkdir -p {params.igcaller}
            wget -O {params.igcaller}/v1.3.zip https://github.com/ferrannadeu/IgCaller/archive/refs/tags/v1.3.zip && unzip {params.igcaller}/v1.3.zip -d {params.igcaller}/ && rm {params.igcaller}/v1.3.zip # MODIFIED to v1.3
        fi

        touch {output.complete}
        """


# Runs IgCaller in paired mode
rule _igcaller_run_paired:
    input:
        installed = CFG["igcaller_path"] + "igcaller-1.3.installed", # MODIFIED to v1.3
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        tsv = CFG["dirs"]["igcaller"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_IgCaller/{tumour_id}_output_filtered.tsv"
    log:
        stdout = CFG["logs"]["igcaller"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/igcaller_run.stdout.log",
        stderr = CFG["logs"]["igcaller"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/igcaller_run.stderr.log"
    params:
        opts = CFG["options"]["_igcaller_run"],
        version = op.switch_on_wildcard("genome_build", CFG["genome_version"]),
        chr_annot = op.switch_on_wildcard("genome_build", CFG["chr_annotation"]),
        seq = op.switch_on_wildcard("seq_type", CFG["switches"]["_igcaller_run"]),
        outdir = CFG["dirs"]["igcaller"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/",
        igcaller = CFG["igcaller_path"]
    wildcard_constraints:
        pair_status = "matched|unmatched"
    conda:
        CFG["conda_envs"]["igcaller"]
    threads:
        CFG["threads"]["_igcaller_run"]
    resources:
        **CFG["resources"]["_igcaller_run"]
    group:
        "input_and_igcaller_run"
    shell:
        op.as_one_line("""
        python3 {params.igcaller}/IgCaller-1.3/IgCaller.py
        -I {params.igcaller}/IgCaller-1.3/IgCaller_reference_files/
        -V {params.version} -T {input.tumour_bam} -N {input.normal_bam} -R {input.fasta}
        -C {params.chr_annot} {params.seq} -o {params.outdir} {params.opts} -@ {threads}
        &&
        rm {params.outdir}/*/*.ba[mi]
        """) # MODIFIED to v1.3 and removed intermediate bam files


# Runs IgCaller in unpaired mode
rule _igcaller_run_unpaired:
    input:
        installed = CFG["igcaller_path"] + "igcaller-1.3.installed", # MODIFIED to v1.3
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        tsv = CFG["dirs"]["igcaller"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_IgCaller/{tumour_id}_output_filtered.tsv"
    log:
        stdout = CFG["logs"]["igcaller"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/igcaller_run.stdout.log",
        stderr = CFG["logs"]["igcaller"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/igcaller_run.stderr.log"
    params:
        opts = CFG["options"]["_igcaller_run"],
        version = op.switch_on_wildcard("genome_build", CFG["genome_version"]),
        chr_annot = op.switch_on_wildcard("genome_build", CFG["chr_annotation"]),
        seq = op.switch_on_wildcard("seq_type", CFG["switches"]["_igcaller_run"]),
        outdir = CFG["dirs"]["igcaller"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/",
        igcaller = CFG["igcaller_path"]
    wildcard_constraints:
        pair_status = "no_normal"
    conda:
        CFG["conda_envs"]["igcaller"]
    threads:
        CFG["threads"]["_igcaller_run"]
    resources:
        **CFG["resources"]["_igcaller_run"]
    group: 
        "input_and_igcaller_run"
    shell:
        op.as_one_line(""" 
        python3 {params.igcaller}/IgCaller-1.3/IgCaller.py 
        -I {params.igcaller}/IgCaller-1.3/IgCaller_reference_files/
        -V {params.version} -T {input.tumour_bam} -R {input.fasta} {params.seq}
        -C {params.chr_annot} -o {params.outdir} {params.opts} -@ {threads}
        &&
        rm {params.outdir}/*/*.ba[mi]
        """) # MODIFIED to v1.3 and removed intermediate bam files


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _igcaller_output_tsv:
    input:
        tsv = CFG["dirs"]["igcaller"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_IgCaller/{tumour_id}_output_filtered.tsv"
    output:
        tsv = CFG["dirs"]["outputs"] + "tsv/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.output.filtered.tsv"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _igcaller_all:
    input:
        expand(
            [
                str(rules._igcaller_output_tsv.output.tsv),
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
