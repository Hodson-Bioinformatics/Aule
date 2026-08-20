#!/usr/bin/env snakemake

import hashlib
import os
import sys

import oncopipe as op

min_oncopipe_version = "1.0.11"
import pkg_resources
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit(
        "The packaging module dependency is missing. Please install it "
        "('pip install packaging') and update oncopipe."
    )

current_version = pkg_resources.get_distribution("oncopipe").version
if version.parse(current_version) < version.parse(min_oncopipe_version):
    sys.exit(
        f"This module requires oncopipe >= {min_oncopipe_version}; "
        f"installed version is {current_version}."
    )

CFG = op.setup_module(
    name = "ichorcna_burden_offtarget",
    version = "1.0",
    subdirectories = [
        "inputs",
        "offtarget_bam",
        "bamCoverage",
        "bigwigToWig",
        "convert_wig",
        "compile_wigs",
        "ichorcna_run",
        "burden",
        "burden_plot",
        "outputs"
    ]
)

localrules:
    _ichorcna_burden_offtarget_input_bam,
    _ichorcna_burden_offtarget_input_pon,
    _ichorcna_burden_offtarget_output,
    _ichorcna_burden_offtarget_aggregate,
    _ichorcna_burden_offtarget_all


md5hash = hashlib.md5()
if workflow.conda_prefix:
    conda_prefix = workflow.conda_prefix
else:
    conda_prefix = os.path.abspath(".snakemake/conda")

md5hash.update(conda_prefix.encode())
with open(CFG["conda_envs"]["ichorcna_run"], "rb") as env_file:
    md5hash.update(env_file.read())
env_hash = md5hash.hexdigest()
ichorcna_env = os.path.join(conda_prefix, f"{env_hash}_")
ichorcna_dir = os.path.join(ichorcna_env, "ichorCNA")
ichorcna_scripts_dir = os.path.join(ichorcna_dir, "scripts")

_get_chrom_wigs_burden = lambda wildcards: expand(
    config["lcr-modules"]["ichorcna_burden_offtarget"]["dirs"]["convert_wig"] +
    "{{seq_type}}--{{genome_build}}/bin{{binSize}}kb/"
    "{{tumour_id}}.bin{{binSize}}kb.{chrom}--fixed.wig",
    chrom = open(reference_files(
        f"genomes/{wildcards.genome_build}/"
        "genome_fasta/main_chromosomes_withY.txt"
    )).read().rstrip("\n").split("\n")
)

_get_burden_input_wig = lambda wildcards: (
    config["lcr-modules"]["ichorcna_burden_offtarget"]["inputs"]["sample_wig"]
    .format(**wildcards)
    if config["lcr-modules"]["ichorcna_burden_offtarget"]["inputs"]
             .get("sample_wig", "")
    else (
        config["lcr-modules"]["ichorcna_burden_offtarget"]["dirs"]["compile_wigs"] +
        "{seq_type}--{genome_build}/bin{binSize}kb/"
        "{tumour_id}.bin{binSize}kb.wig"
    ).format(**wildcards)
)

_ICHOR_OUTPUT_LAYOUT = {
    "corrDepth": ("corrDepth", ".corrDepth.txt"),
    "cna": ("binCNA", ".cna.seg"),
    "segTxt": ("seg_txt", ".seg.txt"),
    "plot": ("plot", "_genomeWide.pdf")
}


def _get_burden_ichor_output(wildcards, output_name):
    output_root = config["lcr-modules"]["ichorcna_burden_offtarget"][
        "inputs"
    ].get("ichorcna_output_dir", "")
    if output_root:
        output_directory, suffix = _ICHOR_OUTPUT_LAYOUT[output_name]
        return os.path.join(
            output_root.format(**wildcards),
            f"{wildcards.seq_type}--{wildcards.genome_build}",
            output_directory,
            f"bin{wildcards.binSize}kb",
            f"{wildcards.tumour_id}{suffix}"
        )
    generated = getattr(
        rules._ichorcna_burden_offtarget_ichorcna_run.output, output_name
    )
    return str(generated).format(**wildcards)


_get_burden_corr_depth = lambda wildcards: _get_burden_ichor_output(
    wildcards, "corrDepth"
)
_get_burden_cna = lambda wildcards: _get_burden_ichor_output(
    wildcards, "cna"
)
_get_burden_seg_txt = lambda wildcards: _get_burden_ichor_output(
    wildcards, "segTxt"
)
_get_burden_plot = lambda wildcards: _get_burden_ichor_output(
    wildcards, "plot"
)

_get_gcwig_burden = lambda wildcards: os.path.join(
    ichorcna_dir,
    "inst",
    "extdata",
    "gc_{}_{}kb.wig".format(
        op.switch_on_wildcard(
            "genome_build",
            config["lcr-modules"]["ichorcna_burden_offtarget"]
                  ["options"]["ichorcna_run"]["genome_ucsc"]
        )(wildcards),
        wildcards.binSize
    )
)

_get_mapwig_burden = lambda wildcards: os.path.join(
    ichorcna_dir,
    "inst",
    "extdata",
    "map_{}_{}kb.wig".format(
        op.switch_on_wildcard(
            "genome_build",
            config["lcr-modules"]["ichorcna_burden_offtarget"]
                  ["options"]["ichorcna_run"]["genome_ucsc"]
        )(wildcards),
        wildcards.binSize
    )
)

_get_centromere_burden = lambda wildcards: os.path.join(
    ichorcna_dir,
    "inst",
    "extdata",
    op.switch_on_wildcard(
        "genome_build",
        config["lcr-modules"]["ichorcna_burden_offtarget"]
              ["options"]["ichorcna_run"]["centromere"]
    )(wildcards)
)

_symlink_burden_outputs = lambda input, output: [
    op.relative_symlink(source, destination, in_module=True)
    for source, destination in zip(input, output)
]


rule _ichorcna_burden_offtarget_install_ichorcna:
    params:
        branch = CFG["options"]["ichorcna_branch"],
        directory = ichorcna_dir,
        lib = lambda wildcards, output: os.path.dirname(output.ichorcna)
    output:
        ichorcna = directory(os.path.join(
            ichorcna_env, "lib", "R", "library", "ichorCNA"
        ))
    log:
        stdout = CFG["logs"]["inputs"] + "install_ichorCNA.stdout.log",
        stderr = CFG["logs"]["inputs"] + "install_ichorCNA.stderr.log"
    conda:
        CFG["conda_envs"]["ichorcna_run"]
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


rule _ichorcna_burden_offtarget_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam.bai"
    log:
        CFG["logs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.log"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)


rule _ichorcna_burden_offtarget_input_pon:
    input:
        pon = CFG["inputs"]["ichorcna_pon"]
    output:
        pon = CFG["dirs"]["inputs"] + "pon/{seq_type}--{genome_build}/bin{binSize}kb/ichorCNA_PON_{genome_build}_bin{binSize}_full_median.rds"
    log:
        CFG["logs"]["inputs"] + "pon/{seq_type}--{genome_build}/bin{binSize}kb.log"
    run:
        op.absolute_symlink(input.pon, output.pon)


rule _ichorcna_burden_offtarget_exclude_panel:
    input:
        bam = str(rules._ichorcna_burden_offtarget_input_bam.output.bam),
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


rule _ichorcna_burden_offtarget_bamCoverage:
    input:
        bam = str(rules._ichorcna_burden_offtarget_exclude_panel.output.bam),
        bai = str(rules._ichorcna_burden_offtarget_exclude_panel.output.bai)
    output:
        bw = CFG["dirs"]["bamCoverage"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb.bw"
    params:
        binSize = "{binSize}",
        opts = CFG["options"]["bamCoverage"]
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
        bamCoverage -b {input.bam} --binSize {params.binSize}000 {params.opts} \
            -o {output.bw} -p {threads}
        """


rule _ichorcna_burden_offtarget_bigwigToWig:
    input:
        bw = str(rules._ichorcna_burden_offtarget_bamCoverage.output.bw)
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
        chrom = ".+(?<!--fixed)"
    log:
        CFG["logs"]["bigwigToWig"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb_{chrom}.log"
    shell:
        "bigWigToWig {input.bw} {output.wig} -chrom={wildcards.chrom}"


rule _ichorcna_burden_offtarget_convert_wig:
    input:
        wig = str(rules._ichorcna_burden_offtarget_bigwigToWig.output.wig)
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
        chrom = ".+(?<!--fixed)"
    log:
        CFG["logs"]["convert_wig"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb.{chrom}--fixed.log"
    shell:
        """
        echo -e "fixedStep chrom={wildcards.chrom} start=1 step={wildcards.binSize}000 span={wildcards.binSize}000 " > {output.wig} &&
        intersectBed -a <(bedops --chop {wildcards.binSize}000 --header {input.wig}) \
            -b {input.wig} -wa -wb | awk '{{print $7}}' >> {output.wig} 2>> {log}
        """


rule _ichorcna_burden_offtarget_compile_wigs:
    input:
        wigs = _get_chrom_wigs_burden
    output:
        wig = CFG["dirs"]["compile_wigs"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb.wig"
    resources:
        **CFG["resources"]["compile_wigs"]
    log:
        CFG["logs"]["compile_wigs"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.log"
    conda:
        CFG["conda_envs"]["convert_wig"]
    shell:
        "cat {input.wigs} > {output.wig} 2> {log}"


rule _ichorcna_burden_offtarget_ichorcna_run:
    input:
        wig = _get_burden_input_wig,
        pon = str(rules._ichorcna_burden_offtarget_input_pon.output.pon),
        ichorCNA = str(rules._ichorcna_burden_offtarget_install_ichorcna.output.ichorcna)
    output:
        corrDepth = CFG["dirs"]["ichorcna_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}.correctedDepth.txt",
        param = CFG["dirs"]["ichorcna_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}.params.txt",
        cna = CFG["dirs"]["ichorcna_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}.cna.seg",
        segTxt = CFG["dirs"]["ichorcna_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}.seg.txt",
        seg = CFG["dirs"]["ichorcna_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}.seg",
        plot = CFG["dirs"]["ichorcna_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}/{tumour_id}_genomeWide.pdf"
    params:
        scripts = ichorcna_scripts_dir,
        outDir = CFG["dirs"]["ichorcna_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}",
        name = "{tumour_id}",
        gcwig = _get_gcwig_burden,
        mapwig = _get_mapwig_burden,
        centromere = _get_centromere_burden,
        exome = CFG["options"]["ichorcna_run"]["exome"],
        genomeStyle = op.switch_on_wildcard(
            "genome_build", CFG["options"]["ichorcna_run"]["genomeStyle"]
        ),
        genomeBuild = "{genome_build}",
        opts = CFG["options"]["ichorcna_run"]["opts"]
    conda:
        CFG["conda_envs"]["ichorcna_run"]
    threads:
        CFG["threads"]["ichorcna_run"]
    resources:
        **CFG["resources"]["ichorcna_run"]
    wildcard_constraints:
        binSize = "|".join(["10", "50", "500", "1000"])
    log:
        stdout = CFG["logs"]["ichorcna_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.stdout.log",
        stderr = CFG["logs"]["ichorcna_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.stderr.log"
    shell:
        op.as_one_line("""
        Rscript {params.scripts}/runIchorCNA.R
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
        {params.opts}
        > {log.stdout} 2> {log.stderr}
        """)


rule _ichorcna_burden_offtarget_burden:
    input:
        bins = _get_burden_corr_depth,
        seg = _get_burden_seg_txt,
        pon = str(rules._ichorcna_burden_offtarget_input_pon.output.pon),
        script = CFG["options"]["src"]["burden_index"]
    output:
        tsv = CFG["dirs"]["burden"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.burden.tsv"
    params:
        sample = "{tumour_id}",
        exclude_bed = CFG["options"]["burden"]["exclude_bed"],
        winsor_p = CFG["options"]["burden"]["winsor_p"],
        window_widths = CFG["options"]["burden"]["window_widths"],
        detect_w = CFG["options"]["burden"]["detect_w"],
        n_perm = CFG["options"]["burden"]["n_perm"],
        n_boot = CFG["options"]["burden"]["n_boot"],
        block = CFG["options"]["burden"]["block"],
        seed = CFG["options"]["burden"]["seed"],
        min_bins = CFG["options"]["burden"]["min_bins"],
        fga_k = CFG["options"]["burden"]["fga_k"],
        sex_flag = "--include_sex_chrom" if CFG["options"]["burden"]["include_sex_chrom"] else ""
    conda:
        CFG["conda_envs"]["ichorcna_run"]
    resources:
        **CFG["resources"]["burden"]
    log:
        CFG["logs"]["burden"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.log"
    shell:
        op.as_one_line("""
        Rscript {input.script}
        --bins {input.bins}
        --seg {input.seg}
        --pon {input.pon}
        --output {output.tsv}
        --sample {params.sample}
        --binsize_kb {wildcards.binSize}
        --exclude_bed '{params.exclude_bed}'
        --min_bins {params.min_bins}
        --winsor_p {params.winsor_p}
        --window_widths {params.window_widths}
        --detect_w {params.detect_w}
        --n_perm {params.n_perm}
        --n_boot {params.n_boot}
        --block {params.block}
        --seed {params.seed}
        --fga_k {params.fga_k}
        {params.sex_flag}
        > {log} 2>&1
        """)


rule _ichorcna_burden_offtarget_plot_vc_diagnostic:
    input:
        table = str(rules._ichorcna_burden_offtarget_burden.output.tsv),
        script = CFG["options"]["src"]["plot_vc_diagnostic"]
    output:
        pdf = CFG["dirs"]["burden_plot"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.vc_diagnostic.pdf"
    params:
        sample = "{tumour_id}"
    conda:
        CFG["conda_envs"]["ichorcna_run"]
    resources:
        **CFG["resources"]["plot_vc_diagnostic"]
    log:
        CFG["logs"]["burden_plot"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.log"
    shell:
        op.as_one_line("""
        Rscript {input.script}
        --table {input.table}
        --sample {params.sample}
        --out {output.pdf}
        > {log} 2>&1
        """)


rule _ichorcna_burden_offtarget_output:
    input:
        corrDepth = _get_burden_corr_depth,
        cna = _get_burden_cna,
        seg_txt = _get_burden_seg_txt,
        plot = _get_burden_plot,
        burden = str(rules._ichorcna_burden_offtarget_burden.output.tsv),
        burden_diagnostic = str(rules._ichorcna_burden_offtarget_plot_vc_diagnostic.output.pdf)
    output:
        corrDepth = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/corrDepth/bin{binSize}kb/{tumour_id}.corrDepth.txt",
        cna = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/binCNA/bin{binSize}kb/{tumour_id}.cna.seg",
        seg_txt = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/seg_txt/bin{binSize}kb/{tumour_id}.seg.txt",
        plot = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/plot/bin{binSize}kb/{tumour_id}_genomeWide.pdf",
        burden = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/burden/bin{binSize}kb/{tumour_id}.burden.tsv",
        burden_diagnostic = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/burden_diagnostic/bin{binSize}kb/{tumour_id}.vc_diagnostic.pdf"
    log:
        CFG["logs"]["outputs"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.log"
    run:
        _symlink_burden_outputs(input, output)


rule _ichorcna_burden_offtarget_aggregate:
    input:
        burdens = expand(
            str(rules._ichorcna_burden_offtarget_output.output.burden),
            zip,
            seq_type = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id = CFG["runs"]["tumour_sample_id"],
            binSize = [CFG["options"]["binsize_kb"]] * len(CFG["runs"]["tumour_sample_id"])
        ),
        script = CFG["options"]["src"]["aggregate_burden"]
    output:
        cohort = CFG["dirs"]["outputs"] + "aggregate/all_samples.burden.cohort_relative.tsv",
        qc = CFG["dirs"]["outputs"] + "aggregate/all_samples.burden.qc.tsv"
    params:
        input_dir = lambda wildcards, output: os.path.dirname(
            os.path.dirname(output.cohort)
        )
    conda:
        CFG["conda_envs"]["ichorcna_run"]
    resources:
        **CFG["resources"]["aggregate"]
    log:
        CFG["logs"]["outputs"] + "aggregate_burden.log"
    shell:
        op.as_one_line("""
        Rscript {input.script}
        --input_dir {params.input_dir}
        --output {output.cohort}
        --qc_output {output.qc}
        > {log} 2>&1
        """)


rule _ichorcna_burden_offtarget_all:
    input:
        expand(
            [
                str(rules._ichorcna_burden_offtarget_output.output.corrDepth),
                str(rules._ichorcna_burden_offtarget_output.output.cna),
                str(rules._ichorcna_burden_offtarget_output.output.seg_txt),
                str(rules._ichorcna_burden_offtarget_output.output.plot),
                str(rules._ichorcna_burden_offtarget_output.output.burden),
                str(rules._ichorcna_burden_offtarget_output.output.burden_diagnostic)
            ],
            zip,
            seq_type = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id = CFG["runs"]["tumour_sample_id"],
            binSize = [CFG["options"]["binsize_kb"]] * len(CFG["runs"]["tumour_sample_id"])
        ) + [
            str(rules._ichorcna_burden_offtarget_aggregate.output.cohort),
            str(rules._ichorcna_burden_offtarget_aggregate.output.qc)
        ]


op.cleanup_module(CFG)
