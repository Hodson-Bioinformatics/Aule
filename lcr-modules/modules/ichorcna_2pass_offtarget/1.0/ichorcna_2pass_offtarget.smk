#!/usr/bin/env snakemake


# ---------------------------------------------------------------------------- #
##### ATTRIBUTION #####
# ---------------------------------------------------------------------------- #

# Module author: Joanna A. Krupka
# Additional contributors: N/A


# ---------------------------------------------------------------------------- #
##### SETUP #####
# ---------------------------------------------------------------------------- #

### Modules ###

import oncopipe as op
import hashlib

# Check that the oncopipe dependency is up-to-date.
min_oncopipe_version = "1.0.11"
import pkg_resources
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it "
             "('pip install packaging') and ensure you are using the most "
             "up-to-date oncopipe version")

current_version = pkg_resources.get_distribution("oncopipe").version
if version.parse(current_version) < version.parse(min_oncopipe_version):
    logger.warning(
        '\x1b[0;31;40m' +
        f'ERROR: oncopipe version installed: {current_version}' + "\n" +
        f'ERROR: This module requires oncopipe version >= {min_oncopipe_version}. '
        'Please update oncopipe in your environment' +
        '\x1b[0m'
    )
    sys.exit("Instructions for updating to the current version of oncopipe are "
             "available at https://lcr-modules.readthedocs.io/en/latest/ "
             "(use option 2)")

### Directories ###
CFG = op.setup_module(
    name           = "ichorcna_2pass_offtarget",
    version        = "1.0",
    subdirectories = [
        "inputs",
        "offtarget_bam",
        "bamCoverage",
        "bigwigToWig",
        "convert_wig",
        "compile_wigs",
        "pass1_run",
        "pass1_qc",
        "pass2_run",
        "confidence",
        "outputs"
    ]
)

localrules:
    _ichorcna_2pass_offtarget_input_bam,
    _ichorcna_2pass_offtarget_input_pon,
    _ichorcna_2pass_offtarget_output,
    _ichorcna_2pass_offtarget_aggregate,
    _ichorcna_2pass_offtarget_all


# ---------------------------------------------------------------------------- #
##### CONDA / ICHORCNA INSTALLATION PATH #####
# ---------------------------------------------------------------------------- #

# Derive the path to the ichorCNA conda environment so the R package
# installation step and script paths can be determined deterministically.
md5hash = hashlib.md5()
if workflow.conda_prefix:
    conda_prefix = workflow.conda_prefix
else:
    conda_prefix = os.path.abspath(".snakemake/conda")

md5hash.update(conda_prefix.encode())
f = open(CFG["conda_envs"]["ichorcna_run"], 'rb')
md5hash.update(f.read())
f.close()
h = md5hash.hexdigest()

ichorcna_env         = conda_prefix + "/" + h + "_"
ichorcna_dir         = ichorcna_env + "/ichorCNA"
ichorcna_scripts_dir = ichorcna_env + "/ichorCNA/scripts/"


# ---------------------------------------------------------------------------- #
##### RULES #####
# ---------------------------------------------------------------------------- #

# ---- Install ichorCNA R package into the conda environment -----------------

rule _ichorcna_2pass_offtarget_install_ichorcna:
    params:
        branch    = CFG['options']['ichorcna_branch'],
        directory = ichorcna_dir,
        lib       = ichorcna_env + "/lib/R/library"
    output:
        ichorcna  = directory(ichorcna_env + "/lib/R/library/ichorCNA")
    log:
        stdout = CFG["logs"]["inputs"] + "install_ichorCNA.stdout.log",
        stderr = CFG["logs"]["inputs"] + "install_ichorCNA.stderr.log"
    conda:
        CFG['conda_envs']['ichorcna_run']
    shell:
        op.as_one_line("""
        rm -rf {params.directory} &&
        git clone --branch {params.branch} --single-branch
            https://github.com/broadinstitute/ichorCNA.git {params.directory} &&
        mkdir -p {params.lib} &&
        R -q -e "options(timeout=9999999);
                 install.packages('{params.directory}',
                 lib='{params.lib}', repos=NULL, type='source')"
        > {log.stdout} 2> {log.stderr}
        """)


# ---- Input symlinks ---------------------------------------------------------

rule _ichorcna_2pass_offtarget_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam.bai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)


rule _ichorcna_2pass_offtarget_input_pon:
    input:
        pon = CFG["inputs"]["ichorcna_pon"]
    output:
        pon = CFG["dirs"]["inputs"] + "pon/{seq_type}--{genome_build}/bin{binSize}kb/ichorCNA_PON_{genome_build}_bin{binSize}_full_median.rds"
    run:
        op.absolute_symlink(input.pon, output.pon)


# ---- Off-target BAM extraction ---------------------------------------------

rule _ichorcna_2pass_offtarget_exclude_panel:
    input:
        bam = str(rules._ichorcna_2pass_offtarget_input_bam.output.bam),
        bed = CFG["inputs"]["panel_bed"]
    output:
        bam = temp(CFG["dirs"]["offtarget_bam"] + "{seq_type}--{genome_build}/{tumour_id}_offtarget.bam"),
        bai = temp(CFG["dirs"]["offtarget_bam"] + "{seq_type}--{genome_build}/{tumour_id}_offtarget.bam.bai")
    conda:
        CFG["conda_envs"]["offtarget_bam"]
    threads:
        CFG["threads"]["offtarget_bam"]
    resources:
        **CFG["resources"]["offtarget_bam"]
    log:
        CFG["logs"]["offtarget_bam"] + "{seq_type}--{genome_build}/{tumour_id}_offtarget.log"
    shell:
        op.as_one_line("""
        bedtools intersect -abam {input.bam} -b {input.bed} -v > {output.bam} 2>> {log} &&
        samtools index {output.bam} 2>> {log}
        """)


# ---- WIG generation (shared by both passes) --------------------------------

rule _ichorcna_2pass_offtarget_bamCoverage:
    input:
        bam = str(rules._ichorcna_2pass_offtarget_exclude_panel.output.bam),
        bai = str(rules._ichorcna_2pass_offtarget_exclude_panel.output.bai)
    output:
        bw = CFG["dirs"]["bamCoverage"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb.bw"
    params:
        binSize = "{binSize}",
        opts    = CFG["options"]["bamCoverage"]
    conda:
        CFG["conda_envs"]["bamCoverage"]
    threads:
        CFG["threads"]["bamCoverage"]
    resources:
        **CFG["resources"]["bamCoverage"]
    wildcard_constraints:
        binSize = "|".join(["10", "50", "500", "1000"])
    log:
        CFG["logs"]["bamCoverage"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb.log"
    shell:
        """
            bamCoverage -b {input.bam} --binSize {params.binSize}000 {params.opts} -o {output.bw} -p {threads}
        """


rule _ichorcna_2pass_offtarget_bigwigToWig:
    input:
        bw = str(rules._ichorcna_2pass_offtarget_bamCoverage.output.bw)
    output:
        wig = temp(CFG["dirs"]["bigwigToWig"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb_{chrom}.wig")
    conda:
        CFG["conda_envs"]["bigwigToWig"]
    threads:
        CFG["threads"]["bigwigToWig"]
    resources:
        **CFG["resources"]["bigwigToWig"]
    wildcard_constraints:
        binSize = "|".join(["10", "50", "500", "1000"]),
        chrom   = ".+(?<!--fixed)"
    log:
        CFG["logs"]["bigwigToWig"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb_{chrom}.log"
    shell:
        """
            bigWigToWig {input.bw} {output.wig} -chrom={wildcards.chrom}
        """


rule _ichorcna_2pass_offtarget_convert_wig:
    input:
        wig = str(rules._ichorcna_2pass_offtarget_bigwigToWig.output.wig)
    output:
        wig = temp(CFG["dirs"]["convert_wig"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb.{chrom}--fixed.wig")
    conda:
        CFG["conda_envs"]["convert_wig"]
    threads:
        CFG["threads"]["convert_wig"]
    resources:
        **CFG["resources"]["convert_wig"]
    wildcard_constraints:
        binSize = "|".join(["10", "50", "500", "1000"]),
        chrom   = ".+(?<!--fixed)"
    log:
        CFG["logs"]["convert_wig"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb.{chrom}--fixed.log"
    shell:
        """
            echo -e "fixedStep chrom={wildcards.chrom} start=1 step={wildcards.binSize}000 span={wildcards.binSize}000 " > {output.wig} &&
            intersectBed -a <( bedops --chop {wildcards.binSize}000 --header  {input.wig} ) -b {input.wig} -wa -wb | awk '{{print $7}}' >> {output.wig} 2>> {log}
        """


def _get_chrom_wigs_2pass(wildcards):
    CFG = config["lcr-modules"]["ichorcna_2pass_offtarget"]
    chrs = reference_files("genomes/" + wildcards.genome_build + "/genome_fasta/main_chromosomes_withY.txt")
    with open(chrs) as fh:
        chrs = fh.read().rstrip("\n").split("\n")
    return expand(
        CFG["dirs"]["convert_wig"] + "{{seq_type}}--{{genome_build}}/bin{{binSize}}kb/{{tumour_id}}.bin{{binSize}}kb.{chrom}--fixed.wig",
        chrom = chrs
    )


rule _ichorcna_2pass_offtarget_compile_wigs:
    input:
        wigs = _get_chrom_wigs_2pass
    output:
        wig = CFG["dirs"]["compile_wigs"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb.wig"
    resources:
        **CFG["resources"]["compile_wigs"]
    shell:
        """
            cat {input.wigs} > {output.wig}
        """


# ---- Optional external WIG input (bypasses BAM -> WIG pipeline) ------------

def _get_2pass_input_wig(wildcards):
    """Return the compiled WIG path.

    If ``sample_wig`` is set in the module inputs (non-empty string), expand
    its wildcards and return it directly — this re-uses a WIG already produced
    by another module (e.g. ichorcna_offtarget) and skips the entire
    BAM-to-WIG pipeline for this module.

    Otherwise fall back to the WIG produced by the compile_wigs rule.
    """
    # Re-bind CFG locally so this function survives op.cleanup_module()
    CFG = config["lcr-modules"]["ichorcna_2pass_offtarget"]
    sample_wig = CFG["inputs"].get("sample_wig", "")
    if sample_wig:
        return sample_wig.format(**wildcards)
    return (
        CFG["dirs"]["compile_wigs"]
        + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb.wig"
    ).format(**wildcards)


# ---- Helper functions for ichorCNA reference files -------------------------

def _get_gcwig_2pass(wildcards):
    genome_ucsc = op.switch_on_wildcard(
        "genome_build",
        config["lcr-modules"]["ichorcna_2pass_offtarget"]["options"]["ichorcna_run"]["genome_ucsc"]
    )(wildcards)
    return ichorcna_dir + "/inst/extdata/gc_" + genome_ucsc + "_" + wildcards.binSize + "kb.wig"

def _get_mapwig_2pass(wildcards):
    genome_ucsc = op.switch_on_wildcard(
        "genome_build",
        config["lcr-modules"]["ichorcna_2pass_offtarget"]["options"]["ichorcna_run"]["genome_ucsc"]
    )(wildcards)
    return ichorcna_dir + "/inst/extdata/map_" + genome_ucsc + "_" + wildcards.binSize + "kb.wig"

def _get_centromere_2pass(wildcards):
    centromere = op.switch_on_wildcard(
        "genome_build",
        config["lcr-modules"]["ichorcna_2pass_offtarget"]["options"]["ichorcna_run"]["centromere"]
    )(wildcards)
    return ichorcna_dir + "/inst/extdata/" + centromere


# ---------------------------------------------------------------------------- #
##### PASS 1 — BROAD EXPLORATION #############################################
# ---------------------------------------------------------------------------- #

rule _ichorcna_2pass_offtarget_pass1_run:
    input:
        wig      = _get_2pass_input_wig,
        pon      = str(rules._ichorcna_2pass_offtarget_input_pon.output.pon),
        ichorCNA = str(rules._ichorcna_2pass_offtarget_install_ichorcna.output.ichorcna)
    output:
        corrDepth = CFG["dirs"]["pass1_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}.correctedDepth.txt",
        param     = CFG["dirs"]["pass1_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}.params.txt",
        cna       = CFG["dirs"]["pass1_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}.cna.seg",
        segTxt    = CFG["dirs"]["pass1_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}.seg.txt",
        seg       = CFG["dirs"]["pass1_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}.seg",
        plot      = CFG["dirs"]["pass1_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}/{tumour_id}_genomeWide.pdf"
    params:
        ichorcna_scripts_dir = ichorcna_scripts_dir,
        outDir               = CFG["dirs"]["pass1_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}",
        name                 = "{tumour_id}",
        gcwig                = _get_gcwig_2pass,
        mapwig               = _get_mapwig_2pass,
        centromere           = _get_centromere_2pass,
        exome                = CFG["options"]["ichorcna_run"]["exome"],
        genomeStyle          = op.switch_on_wildcard("genome_build", CFG["options"]["ichorcna_run"]["genomeStyle"]),
        genomeBuild          = "{genome_build}",
        pass1_opts           = CFG["options"]["ichorcna_run"]["pass1_opts"]
    conda:
        CFG["conda_envs"]["ichorcna_run"]
    threads:
        CFG["threads"]["ichorcna_run"]
    resources:
        **CFG["resources"]["ichorcna_run"]
    wildcard_constraints:
        binSize = "|".join(["10", "50", "500", "1000"])
    log:
        stdout = CFG["logs"]["pass1_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.stdout.log",
        stderr = CFG["logs"]["pass1_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.stderr.log"
    shell:
        op.as_one_line("""
            Rscript {params.ichorcna_scripts_dir}/runIchorCNA.R
            --id {params.name}
            --WIG {input.wig}
            --gcWig {params.gcwig}
            --mapWig {params.mapwig}
            --normalPanel {input.pon}
            --exons.bed {params.exome}
            --genomeStyle {params.genomeStyle}
            --genomeBuild {params.genomeBuild}
            --centromere {params.centromere}
            --outDir {params.outDir}
            {params.pass1_opts}
            > {log.stdout} 2> {log.stderr}
        """)


# ---------------------------------------------------------------------------- #
##### PASS 1 QC — DECISION LOGIC #############################################
# ---------------------------------------------------------------------------- #

rule _ichorcna_2pass_offtarget_pass1_qc:
    input:
        param    = str(rules._ichorcna_2pass_offtarget_pass1_run.output.param),
        qc_script = CFG["options"]["src"]["pass1_qc"]
    output:
        strategy = CFG["dirs"]["pass1_qc"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.strategy.tsv"
    params:
        sample           = "{tumour_id}",
        max_normal_init  = CFG["options"]["ichorcna_run"]["pass1_max_normal_init"]
    conda:
        CFG["conda_envs"]["ichorcna_run"]
    threads: 1
    resources:
        **CFG["resources"]["pass1_qc"]
    log:
        stdout = CFG["logs"]["pass1_qc"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.stdout.log",
        stderr = CFG["logs"]["pass1_qc"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.stderr.log"
    shell:
        op.as_one_line("""
            Rscript {input.qc_script}
            --params  {input.param}
            --strategy {output.strategy}
            --sample  {params.sample}
            --max_normal_init {params.max_normal_init}
            > {log.stdout} 2> {log.stderr}
        """)


# ---------------------------------------------------------------------------- #
##### PASS 2 — STRATEGY-DRIVEN REFINEMENT ####################################
# ---------------------------------------------------------------------------- #

rule _ichorcna_2pass_offtarget_pass2_run:
    input:
        wig              = _get_2pass_input_wig,
        pon              = str(rules._ichorcna_2pass_offtarget_input_pon.output.pon),
        ichorCNA         = str(rules._ichorcna_2pass_offtarget_install_ichorcna.output.ichorcna),
        strategy         = str(rules._ichorcna_2pass_offtarget_pass1_qc.output.strategy),
        pass1_corrDepth  = str(rules._ichorcna_2pass_offtarget_pass1_run.output.corrDepth),
        pass1_param      = str(rules._ichorcna_2pass_offtarget_pass1_run.output.param),
        pass1_cna        = str(rules._ichorcna_2pass_offtarget_pass1_run.output.cna),
        pass1_segTxt     = str(rules._ichorcna_2pass_offtarget_pass1_run.output.segTxt),
        pass1_seg        = str(rules._ichorcna_2pass_offtarget_pass1_run.output.seg),
        pass1_plot       = str(rules._ichorcna_2pass_offtarget_pass1_run.output.plot)
    output:
        corrDepth = CFG["dirs"]["pass2_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}.correctedDepth.txt",
        param     = CFG["dirs"]["pass2_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}.params.txt",
        cna       = CFG["dirs"]["pass2_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}.cna.seg",
        segTxt    = CFG["dirs"]["pass2_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}.seg.txt",
        seg       = CFG["dirs"]["pass2_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}.seg",
        plot      = CFG["dirs"]["pass2_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}/{tumour_id}_genomeWide.pdf"
    params:
        ichorcna_scripts_dir = ichorcna_scripts_dir,
        outDir               = CFG["dirs"]["pass2_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}",
        name                 = "{tumour_id}",
        gcwig                = _get_gcwig_2pass,
        mapwig               = _get_mapwig_2pass,
        centromere           = _get_centromere_2pass,
        exome                = CFG["options"]["ichorcna_run"]["exome"],
        genomeStyle          = op.switch_on_wildcard("genome_build", CFG["options"]["ichorcna_run"]["genomeStyle"]),
        genomeBuild          = "{genome_build}"
    conda:
        CFG["conda_envs"]["ichorcna_run"]
    threads:
        CFG["threads"]["ichorcna_run"]
    resources:
        **CFG["resources"]["ichorcna_run"]
    wildcard_constraints:
        binSize = "|".join(["10", "50", "500", "1000"])
    log:
        stdout = CFG["logs"]["pass2_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.stdout.log",
        stderr = CFG["logs"]["pass2_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.stderr.log"
    shell:
        # Read the strategy from the TSV, then either:
        #   ACCEPT / UNDETECTABLE  → forward Pass 1 outputs unchanged
        #   UNRESOLVABLE           → re-run with strict subclonal suppression
        #   SMOOTH / RESTRICT / PLOIDY_LOCK → re-run with strategy-specific params
        op.as_one_line("""
        STRATEGY=$(grep -P "^strategy\t" {input.strategy} | cut -f2) ;
        if [ "$STRATEGY" = "ACCEPT" ] || [ "$STRATEGY" = "UNDETECTABLE" ] ; then
            cp {input.pass1_corrDepth} {output.corrDepth} &&
            cp {input.pass1_param}     {output.param}     &&
            cp {input.pass1_cna}       {output.cna}       &&
            cp {input.pass1_segTxt}    {output.segTxt}    &&
            cp {input.pass1_seg}       {output.seg}       &&
            mkdir -p $(dirname {output.plot})              &&
            cp {input.pass1_plot}      {output.plot} ;
        elif [ "$STRATEGY" = "UNRESOLVABLE" ] ; then
            Rscript {params.ichorcna_scripts_dir}/runIchorCNA.R
                --id {params.name}
                --WIG {input.wig}
                --gcWig {params.gcwig}
                --mapWig {params.mapwig}
                --normalPanel {input.pon}
                --exons.bed {params.exome}
                --genomeStyle {params.genomeStyle}
                --genomeBuild {params.genomeBuild}
                --centromere {params.centromere}
                --outDir {params.outDir}
                --estimateScPrevalence FALSE
                --scStates "c()"
                --maxFracCNASubclone 0.0
                --maxFracGenomeSubclone 0.0
                --txnE 0.9999
                --txnStrength 10000
                --maxCN 4
                --ploidy "c(2)"
                --lambda "c(1500,1500,1500,1500)"
                --normal "c(0.5,0.6,0.7,0.8,0.9,0.95)"
                > {log.stdout} 2> {log.stderr} ;
        else
            PASS2_NORMAL=$(grep -P "^pass2_normal\t"      {input.strategy} | cut -f2) ;
            PASS2_PLOIDY=$(grep -P "^pass2_ploidy\t"      {input.strategy} | cut -f2) ;
            PASS2_MAXCN=$(grep  -P "^pass2_maxCN\t"       {input.strategy} | cut -f2) ;
            PASS2_TXNE=$(grep   -P "^pass2_txnE\t"        {input.strategy} | cut -f2) ;
            PASS2_TXNSTR=$(grep -P "^pass2_txnStrength\t" {input.strategy} | cut -f2) ;
            Rscript {params.ichorcna_scripts_dir}/runIchorCNA.R
                --id {params.name}
                --WIG {input.wig}
                --gcWig {params.gcwig}
                --mapWig {params.mapwig}
                --normalPanel {input.pon}
                --exons.bed {params.exome}
                --genomeStyle {params.genomeStyle}
                --genomeBuild {params.genomeBuild}
                --centromere {params.centromere}
                --outDir {params.outDir}
                --estimateScPrevalence TRUE
                --scStates "c(1,3)"
                --normal "$PASS2_NORMAL"
                --ploidy "$PASS2_PLOIDY"
                --maxCN $PASS2_MAXCN
                --txnE $PASS2_TXNE
                --txnStrength $PASS2_TXNSTR
                > {log.stdout} 2> {log.stderr} ;
        fi
        """)


# ---------------------------------------------------------------------------- #
##### FINAL CONFIDENCE REPORT ################################################
# ---------------------------------------------------------------------------- #

rule _ichorcna_2pass_offtarget_confidence:
    input:
        pass1_param       = str(rules._ichorcna_2pass_offtarget_pass1_run.output.param),
        pass2_param       = str(rules._ichorcna_2pass_offtarget_pass2_run.output.param),
        strategy          = str(rules._ichorcna_2pass_offtarget_pass1_qc.output.strategy),
        confidence_script = CFG["options"]["src"]["final_confidence"]
    output:
        report = CFG["dirs"]["confidence"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.confidence.tsv"
    params:
        sample = "{tumour_id}"
    conda:
        CFG["conda_envs"]["ichorcna_run"]
    threads: 1
    resources:
        **CFG["resources"]["confidence"]
    log:
        stdout = CFG["logs"]["confidence"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.stdout.log",
        stderr = CFG["logs"]["confidence"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.stderr.log"
    shell:
        op.as_one_line("""
            Rscript {input.confidence_script}
            --pass1_params {input.pass1_param}
            --pass2_params {input.pass2_param}
            --strategy     {input.strategy}
            --output       {output.report}
            --sample       {params.sample}
            > {log.stdout} 2> {log.stderr}
        """)


# ---------------------------------------------------------------------------- #
##### OUTPUTS ################################################################
# ---------------------------------------------------------------------------- #

rule _ichorcna_2pass_offtarget_output:
    input:
        # Pass 2 CNA outputs (Pass 1 is copied here for ACCEPT/UNDETECTABLE)
        corrDepth  = str(rules._ichorcna_2pass_offtarget_pass2_run.output.corrDepth),
        param      = str(rules._ichorcna_2pass_offtarget_pass2_run.output.param),
        cna        = str(rules._ichorcna_2pass_offtarget_pass2_run.output.cna),
        seg_txt    = str(rules._ichorcna_2pass_offtarget_pass2_run.output.segTxt),
        seg        = str(rules._ichorcna_2pass_offtarget_pass2_run.output.seg),
        plot       = str(rules._ichorcna_2pass_offtarget_pass2_run.output.plot),
        # Pass 1 outputs (always retained for QC)
        p1_param   = str(rules._ichorcna_2pass_offtarget_pass1_run.output.param),
        p1_plot    = str(rules._ichorcna_2pass_offtarget_pass1_run.output.plot),
        strategy   = str(rules._ichorcna_2pass_offtarget_pass1_qc.output.strategy),
        confidence = str(rules._ichorcna_2pass_offtarget_confidence.output.report)
    output:
        corrDepth  = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/corrDepth/bin{binSize}kb/{tumour_id}.corrDepth.txt",
        param      = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/param/bin{binSize}kb/{tumour_id}.param.txt",
        cna        = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/binCNA/bin{binSize}kb/{tumour_id}.cna.seg",
        seg_txt    = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/seg_txt/bin{binSize}kb/{tumour_id}.seg.txt",
        seg        = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/seg/bin{binSize}kb/{tumour_id}.seg",
        plot       = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/plot/bin{binSize}kb/{tumour_id}_genomeWide.pdf",
        p1_param   = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/pass1_param/bin{binSize}kb/{tumour_id}.pass1.param.txt",
        p1_plot    = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/pass1_plot/bin{binSize}kb/{tumour_id}_pass1_genomeWide.pdf",
        strategy   = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/strategy/bin{binSize}kb/{tumour_id}.strategy.tsv",
        confidence = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/confidence/bin{binSize}kb/{tumour_id}.confidence.tsv"
    run:
        op.relative_symlink(input.corrDepth,  output.corrDepth,  in_module=True)
        op.relative_symlink(input.param,      output.param,      in_module=True)
        op.relative_symlink(input.cna,        output.cna,        in_module=True)
        op.relative_symlink(input.seg_txt,    output.seg_txt,    in_module=True)
        op.relative_symlink(input.seg,        output.seg,        in_module=True)
        op.relative_symlink(input.plot,       output.plot,       in_module=True)
        op.relative_symlink(input.p1_param,   output.p1_param,   in_module=True)
        op.relative_symlink(input.p1_plot,    output.p1_plot,    in_module=True)
        op.relative_symlink(input.strategy,   output.strategy,   in_module=True)
        op.relative_symlink(input.confidence, output.confidence, in_module=True)


# ---- Aggregation: summary table + plot galleries --------------------------

rule _ichorcna_2pass_offtarget_aggregate:
    input:
        confidence = expand(
            str(rules._ichorcna_2pass_offtarget_output.output.confidence),
            zip,
            seq_type     = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id    = CFG["runs"]["tumour_sample_id"],
            binSize      = [CFG["options"]["binsize_kb"]] * len(CFG["runs"]["tumour_sample_id"])
        ),
        plots = expand(
            str(rules._ichorcna_2pass_offtarget_output.output.plot),
            zip,
            seq_type     = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id    = CFG["runs"]["tumour_sample_id"],
            binSize      = [CFG["options"]["binsize_kb"]] * len(CFG["runs"]["tumour_sample_id"])
        ),
        p1_plots = expand(
            str(rules._ichorcna_2pass_offtarget_output.output.p1_plot),
            zip,
            seq_type     = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id    = CFG["runs"]["tumour_sample_id"],
            binSize      = [CFG["options"]["binsize_kb"]] * len(CFG["runs"]["tumour_sample_id"])
        )
    output:
        summary        = CFG["dirs"]["outputs"] + "aggregate/all_samples.confidence.tsv",
        plots_stamp    = CFG["dirs"]["outputs"] + "aggregate/plots/.done",
        p1_plots_stamp = CFG["dirs"]["outputs"] + "aggregate/pass1_plots/.done"
    run:
        import os
        # ---- Merge all per-sample confidence TSVs into a single table --------
        os.makedirs(os.path.dirname(output.summary), exist_ok=True)
        with open(output.summary, "w") as out_fh:
            header_written = False
            for f in input.confidence:
                with open(f) as fh:
                    lines = fh.readlines()
                    if not header_written:
                        out_fh.writelines(lines)
                        header_written = True
                    else:
                        out_fh.writelines(lines[1:])  # skip header for subsequent files
        # ---- Symlink Pass 2 (final) plots into a flat gallery ---------------
        plots_dir = os.path.dirname(output.plots_stamp)
        os.makedirs(plots_dir, exist_ok=True)
        for src in input.plots:
            dest = os.path.join(plots_dir, os.path.basename(src))
            if os.path.lexists(dest):
                os.remove(dest)
            os.symlink(os.path.realpath(src), dest)
        open(output.plots_stamp, "w").close()
        # ---- Symlink Pass 1 plots into a flat gallery -----------------------
        p1_plots_dir = os.path.dirname(output.p1_plots_stamp)
        os.makedirs(p1_plots_dir, exist_ok=True)
        for src in input.p1_plots:
            dest = os.path.join(p1_plots_dir, os.path.basename(src))
            if os.path.lexists(dest):
                os.remove(dest)
            os.symlink(os.path.realpath(src), dest)
        open(output.p1_plots_stamp, "w").close()


rule _ichorcna_2pass_offtarget_aggregate_params:
    input:
        p1_params = expand(
            str(rules._ichorcna_2pass_offtarget_output.output.p1_param),
            zip,
            seq_type     = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id    = CFG["runs"]["tumour_sample_id"],
            binSize      = [CFG["options"]["binsize_kb"]] * len(CFG["runs"]["tumour_sample_id"])
        ),
        params = expand(
            str(rules._ichorcna_2pass_offtarget_output.output.param),
            zip,
            seq_type     = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id    = CFG["runs"]["tumour_sample_id"],
            binSize      = [CFG["options"]["binsize_kb"]] * len(CFG["runs"]["tumour_sample_id"])
        )
    output:
        tsv = CFG["dirs"]["outputs"] + "aggregate/all_samples_all_models.tsv"
    params:
        r_script = CFG["options"]["src"]["aggregate_params"]
    conda:
        CFG["conda_envs"]["ichorcna_run"]
    threads: 1
    resources:
        **CFG["resources"]["aggregate_params"]
    log:
        stdout = CFG["logs"]["outputs"] + "aggregate_params.stdout.log",
        stderr = CFG["logs"]["outputs"] + "aggregate_params.stderr.log"
    shell:
        op.as_one_line("""
            Rscript {params.r_script}
            {output.tsv}
            {input.p1_params}
            {input.params}
            > {log.stdout} 2> {log.stderr}
        """)


rule _ichorcna_2pass_offtarget_aggregate_pdf:
    input:
        p1_plots = expand(
            str(rules._ichorcna_2pass_offtarget_output.output.p1_plot),
            zip,
            seq_type     = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id    = CFG["runs"]["tumour_sample_id"],
            binSize      = [CFG["options"]["binsize_kb"]] * len(CFG["runs"]["tumour_sample_id"])
        ),
        plots = expand(
            str(rules._ichorcna_2pass_offtarget_output.output.plot),
            zip,
            seq_type     = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id    = CFG["runs"]["tumour_sample_id"],
            binSize      = [CFG["options"]["binsize_kb"]] * len(CFG["runs"]["tumour_sample_id"])
        )
    output:
        pdf = CFG["dirs"]["outputs"] + "aggregate/all_samples_pass1_pass2_genomeWide.pdf"
    params:
        r_script = CFG["options"]["src"]["combine_plots"]
    conda:
        CFG["conda_envs"]["aggregate_plots"]
    threads: 1
    resources:
        **CFG["resources"]["aggregate_pdf"]
    log:
        stdout = CFG["logs"]["outputs"] + "aggregate_pdf.stdout.log",
        stderr = CFG["logs"]["outputs"] + "aggregate_pdf.stderr.log"
    shell:
        op.as_one_line("""
            Rscript {params.r_script}
            {output.pdf}
            {input.p1_plots}
            {input.plots}
            > {log.stdout} 2> {log.stderr}
        """)


rule _ichorcna_2pass_offtarget_all:
    input:
        expand(
            [
                str(rules._ichorcna_2pass_offtarget_output.output.corrDepth),
                str(rules._ichorcna_2pass_offtarget_output.output.param),
                str(rules._ichorcna_2pass_offtarget_output.output.cna),
                str(rules._ichorcna_2pass_offtarget_output.output.seg_txt),
                str(rules._ichorcna_2pass_offtarget_output.output.seg),
                str(rules._ichorcna_2pass_offtarget_output.output.plot),
                str(rules._ichorcna_2pass_offtarget_output.output.p1_param),
                str(rules._ichorcna_2pass_offtarget_output.output.p1_plot),
                str(rules._ichorcna_2pass_offtarget_output.output.strategy),
                str(rules._ichorcna_2pass_offtarget_output.output.confidence)
            ],
            zip,
            seq_type    = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id   = CFG["runs"]["tumour_sample_id"],
            binSize     = [CFG["options"]["binsize_kb"]] * len(CFG["runs"]["tumour_sample_id"])
        ) + [
            str(rules._ichorcna_2pass_offtarget_aggregate.output.summary),
            str(rules._ichorcna_2pass_offtarget_aggregate.output.plots_stamp),
            str(rules._ichorcna_2pass_offtarget_aggregate.output.p1_plots_stamp),
            str(rules._ichorcna_2pass_offtarget_aggregate_params.output.tsv),
            str(rules._ichorcna_2pass_offtarget_aggregate_pdf.output.pdf)
        ]


##### CLEANUP #####

op.cleanup_module(CFG)
