#!/usr/bin/env snakemake

import oncopipe as op


CFG = op.setup_module(
    name = "liquid_biopsy_qc",
    version = "1.0",
    subdirectories = ["inputs", "bam_metrics", "coverage", "per_sample", "qc", "outputs"]
)

SKIP_BAM = bool(CFG["options"].get("skip_bam_metrics", False))
ASSAY_MODE = str(CFG["options"].get("assay_mode", "capture")).lower()
COMPUTE_COVERAGE = bool(CFG["options"].get("compute_coverage", True))
IS_CAPTURE = ASSAY_MODE == "capture"
BIN_SIZE = CFG["options"]["binsize_kb"]

if ASSAY_MODE not in {"capture", "swgs"}:
    raise ValueError("liquid_biopsy_qc options.assay_mode must be capture or swgs")
if IS_CAPTURE and not CFG["options"].get("panel_bed"):
    raise ValueError("liquid_biopsy_qc capture mode requires options.panel_bed")
if not SKIP_BAM and not COMPUTE_COVERAGE:
    required_sidecars = [
        "offtarget_stats", "offtarget_flagstat",
        "offtarget_mosdepth_summary", "offtarget_mosdepth_regions"
    ]
    missing_sidecars = [key for key in required_sidecars
                        if not CFG["inputs"].get(key)]
    if missing_sidecars:
        raise ValueError(
            "liquid_biopsy_qc compute_coverage=False requires: " +
            ", ".join(missing_sidecars)
        )


localrules:
    _liquid_biopsy_qc_pon_reference,
    _liquid_biopsy_qc_all


def _format_input(key, wildcards):
    module_cfg = config["lcr-modules"]["liquid_biopsy_qc"]
    return module_cfg["inputs"][key].format(
        seq_type=wildcards.seq_type,
        genome_build=wildcards.genome_build,
        tumour_id=wildcards.tumour_id,
        sample_id=wildcards.tumour_id,
        binSize=wildcards.binSize
    )


def _panel_dependency(wildcards):
    module_cfg = config["lcr-modules"]["liquid_biopsy_qc"]
    mode = str(module_cfg["options"].get("assay_mode", "capture")).lower()
    return [module_cfg["options"]["panel_bed"]] if mode == "capture" else []


def _bam_dependencies(wildcards):
    if SKIP_BAM:
        return []
    module_cfg = config["lcr-modules"]["liquid_biopsy_qc"]
    prefix = (module_cfg["dirs"]["bam_metrics"] +
              "{seq_type}--{genome_build}/{tumour_id}").format(**wildcards)
    return [
        prefix + suffix for suffix in [
            ".samtools_stats.txt", ".markduplicates.txt",
            ".insert_size_metrics.txt", ".insert_size_histogram.pdf",
            ".hs_metrics.txt", ".passing_reads.txt"
        ]
    ]


COVERAGE_SUFFIXES = {
    "stats": ".coverage.samtools_stats.txt",
    "flagstat": ".coverage.flagstat.txt",
    "summary": ".coverage.mosdepth.summary.txt",
    "regions": ".coverage.regions.bed.gz",
}
EXTERNAL_COVERAGE_KEYS = {
    "stats": "offtarget_stats",
    "flagstat": "offtarget_flagstat",
    "summary": "offtarget_mosdepth_summary",
    "regions": "offtarget_mosdepth_regions",
}


def _coverage_prefix(wildcards):
    module_cfg = config["lcr-modules"]["liquid_biopsy_qc"]
    return (module_cfg["dirs"]["coverage"] +
            "{seq_type}--{genome_build}/{tumour_id}").format(**wildcards)


def _coverage_path(kind, wildcards):
    if COMPUTE_COVERAGE:
        return _coverage_prefix(wildcards) + COVERAGE_SUFFIXES[kind]
    return _format_input(EXTERNAL_COVERAGE_KEYS[kind], wildcards)


def _coverage_dependencies(wildcards):
    if SKIP_BAM:
        return []
    return [_coverage_path(kind, wildcards)
            for kind in ["stats", "flagstat", "summary", "regions"]]


def _coverage_pattern(kind):
    if COMPUTE_COVERAGE:
        return (CFG["dirs"]["coverage"] +
                "{seq_type}--{genome_build}/{tumour_id}" + COVERAGE_SUFFIXES[kind])
    return CFG["inputs"][EXTERNAL_COVERAGE_KEYS[kind]]


def _optional_metadata(wildcards):
    path = config["lcr-modules"]["liquid_biopsy_qc"]["inputs"].get(
        "sample_metadata", ""
    )
    return [path] if path else []


def _coverage_arg(kind, wildcards):
    return "" if SKIP_BAM else _coverage_path(kind, wildcards)


def _provenance_args():
    provenance = CFG["options"]["provenance"]
    return ";".join("{}={}".format(key, value)
                    for key, value in provenance.items())


def _threshold_args():
    thresholds = CFG["options"]["thresholds"]
    return ";".join("{}={}".format(key, value)
                    for key, value in thresholds.items())


def _hs_command(wildcards, input, output, resources):
    java_heap_mb = max(512, int(resources.mem_mb * 0.75))
    if IS_CAPTURE:
        return (
            "picard -Xmx{mem}m BedToIntervalList I={bed} O={interval} SD={dictionary} && "
            "picard -Xmx{mem}m CollectHsMetrics I={bam} O={hs} R={fasta} "
            "TI={interval} BI={interval} VALIDATION_STRINGENCY=LENIENT "
            "LEVEL=ALL_READS CLIP_OVERLAPPING_READS=false"
        ).format(mem=java_heap_mb, bed=input.bed, interval=output.interval_list,
                 dictionary=input.dictionary, bam=input.bam, hs=output.hs,
                 fasta=input.fasta)
    return (
        "echo not_applicable_swgs > {hs} && "
        "echo not_applicable_swgs > {interval}"
    ).format(hs=output.hs, interval=output.interval_list)


if not SKIP_BAM and COMPUTE_COVERAGE:
    if IS_CAPTURE:
        rule _liquid_biopsy_qc_coverage:
            group:
                "liquid_biopsy_qc_sample"
            input:
                bam = CFG["inputs"]["sample_bam"],
                bai = CFG["inputs"]["sample_bai"],
                bed = CFG["options"]["panel_bed"]
            output:
                bam = temp(CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{tumour_id}.offtarget.bam"),
                bai = temp(CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{tumour_id}.offtarget.bam.bai"),
                stats = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{tumour_id}.coverage.samtools_stats.txt",
                flagstat = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{tumour_id}.coverage.flagstat.txt",
                summary = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{tumour_id}.coverage.mosdepth.summary.txt",
                regions = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{tumour_id}.coverage.regions.bed.gz"
            params:
                prefix = lambda w: _coverage_prefix(w) + ".coverage",
                window_bp = CFG["options"]["coverage_bin_size_bp"],
                exclude_flags = CFG["options"]["coverage_exclude_flags"]
            conda:
                CFG["conda_envs"]["coverage"]
            threads:
                CFG["threads"]["coverage"]
            resources:
                **CFG["resources"]["coverage"]
            log:
                CFG["logs"]["coverage"] + "{seq_type}--{genome_build}/{tumour_id}.log"
            shell:
                op.as_one_line("""
                bedtools intersect -abam {input.bam} -b {input.bed} -v > {output.bam} 2> {log} &&
                samtools index -@ {threads} {output.bam} {output.bai} 2>> {log} &&
                samtools stats {output.bam} > {output.stats} 2>> {log} &&
                samtools flagstat {output.bam} > {output.flagstat} 2>> {log} &&
                mosdepth --threads {threads} --flag {params.exclude_flags}
                --by {params.window_bp} --no-per-base {params.prefix} {output.bam} 2>> {log}
                """)
    else:
        rule _liquid_biopsy_qc_coverage:
            group:
                "liquid_biopsy_qc_sample"
            input:
                bam = CFG["inputs"]["sample_bam"],
                bai = CFG["inputs"]["sample_bai"]
            output:
                stats = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{tumour_id}.coverage.samtools_stats.txt",
                flagstat = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{tumour_id}.coverage.flagstat.txt",
                summary = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{tumour_id}.coverage.mosdepth.summary.txt",
                regions = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{tumour_id}.coverage.regions.bed.gz"
            params:
                prefix = lambda w: _coverage_prefix(w) + ".coverage",
                window_bp = CFG["options"]["coverage_bin_size_bp"],
                exclude_flags = CFG["options"]["coverage_exclude_flags"]
            conda:
                CFG["conda_envs"]["coverage"]
            threads:
                CFG["threads"]["coverage"]
            resources:
                **CFG["resources"]["coverage"]
            log:
                CFG["logs"]["coverage"] + "{seq_type}--{genome_build}/{tumour_id}.log"
            shell:
                op.as_one_line("""
                samtools stats {input.bam} > {output.stats} 2> {log} &&
                samtools flagstat {input.bam} > {output.flagstat} 2>> {log} &&
                mosdepth --threads {threads} --flag {params.exclude_flags}
                --by {params.window_bp} --no-per-base {params.prefix} {input.bam} 2>> {log}
                """)


if not SKIP_BAM:
    rule _liquid_biopsy_qc_bam_metrics:
        group:
            "liquid_biopsy_qc_sample"
        input:
            bam = CFG["inputs"]["sample_bam"],
            bai = CFG["inputs"]["sample_bai"],
            bed = _panel_dependency,
            fasta = CFG["inputs"]["reference_fasta"],
            dictionary = CFG["inputs"]["reference_dict"],
            coverage = _coverage_dependencies
        output:
            stats = CFG["dirs"]["bam_metrics"] + "{seq_type}--{genome_build}/{tumour_id}.samtools_stats.txt",
            markdup = CFG["dirs"]["bam_metrics"] + "{seq_type}--{genome_build}/{tumour_id}.markduplicates.txt",
            insert_metrics = CFG["dirs"]["bam_metrics"] + "{seq_type}--{genome_build}/{tumour_id}.insert_size_metrics.txt",
            histogram = CFG["dirs"]["bam_metrics"] + "{seq_type}--{genome_build}/{tumour_id}.insert_size_histogram.pdf",
            hs = CFG["dirs"]["bam_metrics"] + "{seq_type}--{genome_build}/{tumour_id}.hs_metrics.txt",
            passing = CFG["dirs"]["bam_metrics"] + "{seq_type}--{genome_build}/{tumour_id}.passing_reads.txt",
            interval_list = temp(CFG["dirs"]["bam_metrics"] + "{seq_type}--{genome_build}/{tumour_id}.panel.interval_list"),
            marked_bam = temp(CFG["dirs"]["bam_metrics"] + "{seq_type}--{genome_build}/{tumour_id}.marked.bam")
        params:
            hs_command = _hs_command,
            java_heap_mb = lambda wildcards, resources: max(
                512, int(resources.mem_mb * 0.75)
            )
        conda:
            CFG["conda_envs"]["bam_metrics"]
        threads:
            CFG["threads"]["bam_metrics"]
        resources:
            **CFG["resources"]["bam_metrics"]
        log:
            CFG["logs"]["bam_metrics"] + "{seq_type}--{genome_build}/{tumour_id}.log"
        shell:
            op.as_one_line("""
            samtools stats {input.bam} > {output.stats} 2> {log} &&
            samtools view -c -q 30 -F 1796 {input.bam} > {output.passing} 2>> {log} &&
            picard -Xmx{params.java_heap_mb}m MarkDuplicates
            I={input.bam} O={output.marked_bam} M={output.markdup}
            VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=false
            SORTING_COLLECTION_SIZE_RATIO=0.15 2>> {log} &&
            picard -Xmx{params.java_heap_mb}m CollectInsertSizeMetrics
            I={input.bam} O={output.insert_metrics} H={output.histogram}
            VALIDATION_STRINGENCY=LENIENT 2>> {log} &&
            {params.hs_command} 2>> {log}
            """)


rule _liquid_biopsy_qc_collect:
    group:
        "liquid_biopsy_qc_sample"
    input:
        corrected_depth = CFG["inputs"]["corrected_depth"],
        params = CFG["inputs"]["params"],
        burden = CFG["inputs"]["burden"],
        pon = CFG["inputs"]["pon"],
        gc_wig = CFG["inputs"]["gc_wig"],
        map_wig = CFG["inputs"]["map_wig"],
        bam_metrics = _bam_dependencies,
        coverage = _coverage_dependencies,
        metadata = _optional_metadata,
        script = CFG["options"]["src"]["collect"]
    output:
        tsv = CFG["dirs"]["per_sample"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.qc.tsv"
    params:
        bam_prefix = lambda w: (
            config["lcr-modules"]["liquid_biopsy_qc"]["dirs"]["bam_metrics"] +
                                "{seq_type}--{genome_build}/{tumour_id}").format(**w),
        coverage_stats = lambda w: _coverage_arg("stats", w),
        coverage_flagstat = lambda w: _coverage_arg("flagstat", w),
        coverage_summary = lambda w: _coverage_arg("summary", w),
        coverage_regions = lambda w: _coverage_arg("regions", w),
        metadata = CFG["inputs"].get("sample_metadata", ""),
        skip_bam = str(SKIP_BAM).lower(),
        assay_mode = ASSAY_MODE,
        normal_grid_max = CFG["options"]["normal_grid_max"],
        tolerance = CFG["options"]["offtarget_fraction_tolerance"],
        provenance = _provenance_args(),
        thresholds = _threshold_args()
    conda:
        CFG["conda_envs"]["r"]
    resources:
        **CFG["resources"]["collect"]
    log:
        CFG["logs"]["per_sample"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.log"
    shell:
        op.as_one_line("""
        Rscript {input.script}
        --sample '{wildcards.tumour_id}'
        --genome-build '{wildcards.genome_build}'
        --corrected-depth {input.corrected_depth}
        --params {input.params}
        --burden {input.burden}
        --pon {input.pon}
        --gc-wig {input.gc_wig}
        --map-wig {input.map_wig}
        --bam-prefix '{params.bam_prefix}'
        --coverage-stats '{params.coverage_stats}'
        --coverage-flagstat '{params.coverage_flagstat}'
        --coverage-summary '{params.coverage_summary}'
        --coverage-regions '{params.coverage_regions}'
        --metadata '{params.metadata}'
        --skip-bam {params.skip_bam}
        --assay-mode {params.assay_mode}
        --normal-grid-max {params.normal_grid_max}
        --offtarget-fraction-tolerance {params.tolerance}
        --provenance '{params.provenance}'
        --thresholds '{params.thresholds}'
        --binsize-kb {wildcards.binSize}
        --output {output.tsv}
        > {log} 2>&1
        """)


rule _liquid_biopsy_qc_aggregate:
    input:
        samples = expand(
            str(rules._liquid_biopsy_qc_collect.output.tsv),
            zip,
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            binSize=[BIN_SIZE] * len(CFG["runs"]["tumour_sample_id"])
        ),
        pon_reference = CFG["inputs"]["pon_noise_reference"],
        script = CFG["options"]["src"]["aggregate"]
    output:
        per_sample = CFG["dirs"]["outputs"] + "qc/per_sample_qc.tsv",
        summary = CFG["dirs"]["outputs"] + "qc/cohort_qc_summary.tsv",
        correlation = CFG["dirs"]["outputs"] + "qc/burden_correlation_check.tsv",
        provenance = CFG["dirs"]["outputs"] + "qc/provenance.tsv",
        insert_pdf = CFG["dirs"]["outputs"] + "qc/figures/insert_size_distributions.pdf",
        insert_png = CFG["dirs"]["outputs"] + "qc/figures/insert_size_distributions.png",
        noise_pdf = CFG["dirs"]["outputs"] + "qc/figures/sigma_noise_vs_pon.pdf",
        noise_png = CFG["dirs"]["outputs"] + "qc/figures/sigma_noise_vs_pon.png",
        gc_pdf = CFG["dirs"]["outputs"] + "qc/figures/resid_gc_cor_vs_burden_amp.pdf",
        gc_png = CFG["dirs"]["outputs"] + "qc/figures/resid_gc_cor_vs_burden_amp.png",
        flags_pdf = CFG["dirs"]["outputs"] + "qc/figures/flag_frequency.pdf",
        flags_png = CFG["dirs"]["outputs"] + "qc/figures/flag_frequency.png"
    params:
        inputs = lambda wildcards, input: ",".join(input.samples),
        raw_dir = CFG["dirs"]["bam_metrics"],
        pon_reference = lambda wildcards, input: input.pon_reference,
        rho = CFG["options"]["burden_correlation_abs_rho"],
        figure_dir = CFG["dirs"]["outputs"] + "qc/figures"
    conda:
        CFG["conda_envs"]["r"]
    resources:
        **CFG["resources"]["aggregate"]
    log:
        CFG["logs"]["qc"] + "aggregate.log"
    shell:
        op.as_one_line("""
        Rscript {input.script}
        --inputs '{params.inputs}' --per-sample {output.per_sample}
        --summary {output.summary} --correlation {output.correlation}
        --provenance {output.provenance} --rho-threshold {params.rho}
        --raw-dir '{params.raw_dir}' --pon-reference {params.pon_reference}
        --figure-dir {params.figure_dir}
        > {log} 2>&1
        """)


rule _liquid_biopsy_qc_pon_reference:
    input:
        reference = CFG["inputs"]["pon_noise_reference"]
    output:
        reference = CFG["dirs"]["outputs"] + "qc/pon_noise_reference.tsv"
    run:
        op.relative_symlink(input.reference, output.reference, in_module=True)


if not SKIP_BAM:
    rule _liquid_biopsy_qc_multiqc:
        input:
            full = expand(
                CFG["dirs"]["bam_metrics"] + "{seq_type}--{genome_build}/{tumour_id}.samtools_stats.txt",
                zip,
                seq_type=CFG["runs"]["tumour_seq_type"],
                genome_build=CFG["runs"]["tumour_genome_build"],
                tumour_id=CFG["runs"]["tumour_sample_id"]
            ),
            coverage_stats = expand(
                _coverage_pattern("stats"), zip,
                seq_type=CFG["runs"]["tumour_seq_type"],
                genome_build=CFG["runs"]["tumour_genome_build"],
                tumour_id=CFG["runs"]["tumour_sample_id"],
                binSize=[BIN_SIZE] * len(CFG["runs"]["tumour_sample_id"])
            ),
            coverage_summary = expand(
                _coverage_pattern("summary"), zip,
                seq_type=CFG["runs"]["tumour_seq_type"],
                genome_build=CFG["runs"]["tumour_genome_build"],
                tumour_id=CFG["runs"]["tumour_sample_id"],
                binSize=[BIN_SIZE] * len(CFG["runs"]["tumour_sample_id"])
            )
        output:
            html = CFG["dirs"]["outputs"] + "qc/multiqc_report.html"
        params:
            raw_dir = CFG["dirs"]["bam_metrics"],
            out_dir = CFG["dirs"]["outputs"] + "qc"
        conda:
            CFG["conda_envs"]["multiqc"]
        threads:
            CFG["threads"]["multiqc"]
        resources:
            **CFG["resources"]["multiqc"]
        log:
            CFG["logs"]["qc"] + "multiqc.log"
        shell:
            "multiqc --force --outdir {params.out_dir} --filename multiqc_report.html {params.raw_dir} {input.coverage_stats} {input.coverage_summary} > {log} 2>&1"


rule _liquid_biopsy_qc_all:
    input:
        str(rules._liquid_biopsy_qc_aggregate.output.per_sample),
        str(rules._liquid_biopsy_qc_aggregate.output.summary),
        str(rules._liquid_biopsy_qc_aggregate.output.correlation),
        str(rules._liquid_biopsy_qc_aggregate.output.provenance),
        str(rules._liquid_biopsy_qc_pon_reference.output.reference),
        *([] if SKIP_BAM else [str(rules._liquid_biopsy_qc_multiqc.output.html)]),
        str(rules._liquid_biopsy_qc_aggregate.output.insert_pdf),
        str(rules._liquid_biopsy_qc_aggregate.output.insert_png),
        str(rules._liquid_biopsy_qc_aggregate.output.noise_pdf),
        str(rules._liquid_biopsy_qc_aggregate.output.noise_png),
        str(rules._liquid_biopsy_qc_aggregate.output.gc_pdf),
        str(rules._liquid_biopsy_qc_aggregate.output.gc_png),
        str(rules._liquid_biopsy_qc_aggregate.output.flags_pdf),
        str(rules._liquid_biopsy_qc_aggregate.output.flags_png)


op.cleanup_module(CFG)
