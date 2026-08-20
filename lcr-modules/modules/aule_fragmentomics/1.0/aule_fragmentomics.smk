#!/usr/bin/env snakemake

import oncopipe as op


CFG = op.setup_module(
    name = "aule_fragmentomics",
    version = "1.0",
    subdirectories = ["inputs", "per_sample", "outputs"]
)

ENABLE_ALLELE = bool(CFG["options"].get("enable_allele_metrics", True))
ENABLE_WEIGHTING = bool(CFG["options"].get("enable_size_weighting", False))
ENABLE_DELFI = bool(CFG["options"].get("enable_delfi", False))
ENABLE_GRIFFIN = bool(CFG["options"].get("enable_griffin", False))
GENERATE_DELFI_BINS = bool(CFG["options"].get("generate_delfi_bins", True))
PREPARE_GRIFFIN_SITES = bool(CFG["options"].get("prepare_griffin_sites", True))
BUILD_DELFI_REFERENCE = bool(CFG["options"].get("build_delfi_reference", False))
BUILD_GRIFFIN_REFERENCE = bool(CFG["options"].get("build_griffin_reference", False))
REFERENCE_SAMPLE_IDS = list(CFG["options"].get("reference_sample_ids", []))
REFERENCE_MIN_SAMPLES = int(CFG["options"].get("reference_min_samples", 10))

DELFI_BINS_GENERATED = CFG["dirs"]["outputs"] + "assets/delfi_bins.bed"
GRIFFIN_SITES_GENERATED = CFG["dirs"]["outputs"] + "assets/griffin_sites.bed"
DELFI_REFERENCE_GENERATED = CFG["dirs"]["outputs"] + "assets/delfi_healthy_reference.tsv"
GRIFFIN_REFERENCE_GENERATED = CFG["dirs"]["outputs"] + "assets/griffin_healthy_reference.tsv"


def _control_inputs(key):
    pattern = CFG["inputs"].get(key, "")
    return ([pattern.format(sample_id=sample_id) for sample_id in REFERENCE_SAMPLE_IDS]
            if pattern else [])


DELFI_REFERENCE_INPUTS = _control_inputs("delfi_reference_inputs")
GRIFFIN_REFERENCE_INPUTS = _control_inputs("griffin_reference_inputs")

if ENABLE_WEIGHTING and not ENABLE_ALLELE:
    raise ValueError("enable_size_weighting requires enable_allele_metrics")
if ENABLE_DELFI and not GENERATE_DELFI_BINS and not CFG["inputs"].get("delfi_bins_bed"):
    raise ValueError("enable_delfi requires generated bins or inputs.delfi_bins_bed")
if ENABLE_GRIFFIN and PREPARE_GRIFFIN_SITES and not CFG["inputs"].get("griffin_raw_sites_bed"):
    raise ValueError("prepared Griffin sites require inputs.griffin_raw_sites_bed")
if ENABLE_GRIFFIN and not PREPARE_GRIFFIN_SITES and not CFG["inputs"].get("griffin_sites_bed"):
    raise ValueError("enable_griffin requires prepared sites or inputs.griffin_sites_bed")
if (BUILD_DELFI_REFERENCE or BUILD_GRIFFIN_REFERENCE) and REFERENCE_MIN_SAMPLES < 2:
    raise ValueError("reference_min_samples must be at least two")
if BUILD_DELFI_REFERENCE and len(DELFI_REFERENCE_INPUTS) < REFERENCE_MIN_SAMPLES:
    raise ValueError("DELFI reference building requires explicit reference_sample_ids "
                     "and delfi_reference_inputs for at least reference_min_samples controls")
if BUILD_GRIFFIN_REFERENCE and len(GRIFFIN_REFERENCE_INPUTS) < REFERENCE_MIN_SAMPLES:
    raise ValueError("Griffin reference building requires explicit reference_sample_ids "
                     "and griffin_reference_inputs for at least reference_min_samples controls")


localrules:
    _aule_fragmentomics_input_bam,
    _aule_fragmentomics_output,
    _aule_fragmentomics_all


def _formatted_input(key, wildcards):
    value = config["lcr-modules"]["aule_fragmentomics"]["inputs"].get(key, "")
    return value.format(**wildcards) if value else ""


def _variant_dependencies(wildcards):
    return [_formatted_input("variant_vcf", wildcards)] if ENABLE_ALLELE else []


def _weight_dependencies(wildcards):
    return [_formatted_input("size_weights", wildcards)] if ENABLE_WEIGHTING else []


def _regional_dependency(key, enabled, wildcards):
    value = _formatted_input(key, wildcards) if enabled else ""
    return [value] if value else []


def _optional_arg(key, enabled, wildcards):
    return _formatted_input(key, wildcards) if enabled else ""


def _griffin_allowlist_dependencies(wildcards):
    value = _formatted_input("griffin_tf_allowlist", wildcards)
    return [value] if value else []


def _griffin_raw_site_dependencies(wildcards):
    value = _formatted_input("griffin_raw_sites_bed", wildcards)
    return [value] if value and ENABLE_GRIFFIN and PREPARE_GRIFFIN_SITES else []


def _delfi_bins(wildcards):
    return (DELFI_BINS_GENERATED if GENERATE_DELFI_BINS
            else _formatted_input("delfi_bins_bed", wildcards))


def _griffin_sites(wildcards):
    return (GRIFFIN_SITES_GENERATED if PREPARE_GRIFFIN_SITES
            else _formatted_input("griffin_sites_bed", wildcards))


def _delfi_reference(wildcards):
    return (DELFI_REFERENCE_GENERATED if BUILD_DELFI_REFERENCE
            else _formatted_input("delfi_reference", wildcards))


def _griffin_reference(wildcards):
    return (GRIFFIN_REFERENCE_GENERATED if BUILD_GRIFFIN_REFERENCE
            else _formatted_input("griffin_reference", wildcards))


def _selected_dependency(selector, enabled, wildcards):
    value = selector(wildcards) if enabled else ""
    return [value] if value else []


rule _aule_fragmentomics_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam.bai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)


rule _aule_fragmentomics_collect:
    input:
        bam = str(rules._aule_fragmentomics_input_bam.output.bam),
        bai = str(rules._aule_fragmentomics_input_bam.output.bai),
        fasta = CFG["inputs"]["reference_fasta"],
        fasta_fai = CFG["inputs"]["reference_fasta_fai"],
        blacklist = CFG["inputs"]["blacklist_bed"],
        vcf = _variant_dependencies,
        weights = _weight_dependencies,
        script = CFG["options"]["src"]["collect"]
    output:
        summary = CFG["dirs"]["per_sample"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.fragmentomics.tsv",
        size_hist = CFG["dirs"]["per_sample"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.fragment_size_hist.tsv",
        motifs = CFG["dirs"]["per_sample"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.end_motifs.tsv",
        allele_hist = CFG["dirs"]["per_sample"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.allele_fragment_size_hist.tsv",
        variants = CFG["dirs"]["per_sample"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.per_variant.tsv"
    params:
        assay_partition = CFG["options"]["assay_partition"],
        vcf = lambda w: _optional_arg("variant_vcf", ENABLE_ALLELE, w),
        weights = lambda w: _optional_arg("size_weights", ENABLE_WEIGHTING, w),
        pass_only = str(CFG["options"]["pass_variants_only"]).lower(),
        min_mapq = CFG["options"]["min_mapq"],
        min_baseq = CFG["options"]["min_baseq"],
        exclude_flags = CFG["options"]["exclude_flags"],
        min_length = CFG["options"]["min_fragment_length"],
        max_length = CFG["options"]["max_fragment_length"],
        short_min = CFG["options"]["short_min"],
        short_max = CFG["options"]["short_max"],
        long_min = CFG["options"]["long_min"],
        long_max = CFG["options"]["long_max"],
        periodicity_max_length = CFG["options"]["periodicity_max_length"],
        periodicity_min_bp = CFG["options"]["periodicity_min_bp"],
        periodicity_max_bp = CFG["options"]["periodicity_max_bp"],
        min_mutant_fragments = CFG["options"]["min_mutant_fragments"]
    conda:
        CFG["conda_envs"]["python"]
    threads:
        CFG["threads"]["collect"]
    resources:
        **CFG["resources"]["collect"]
    log:
        CFG["logs"]["per_sample"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.log"
    shell:
        op.as_one_line("""
        python {input.script}
        --sample-id '{wildcards.tumour_id}'
        --assay-partition '{params.assay_partition}'
        --bam {input.bam}
        --fasta {input.fasta}
        --blacklist {input.blacklist}
        --variant-vcf '{params.vcf}'
        --size-weights '{params.weights}'
        --pass-variants-only {params.pass_only}
        --min-mapq {params.min_mapq}
        --min-baseq {params.min_baseq}
        --exclude-flags {params.exclude_flags}
        --min-length {params.min_length}
        --max-length {params.max_length}
        --short-min {params.short_min}
        --short-max {params.short_max}
        --long-min {params.long_min}
        --long-max {params.long_max}
        --periodicity-max-length {params.periodicity_max_length}
        --periodicity-min-bp {params.periodicity_min_bp}
        --periodicity-max-bp {params.periodicity_max_bp}
        --min-mutant-fragments {params.min_mutant_fragments}
        --summary {output.summary}
        --size-hist {output.size_hist}
        --motifs {output.motifs}
        --allele-hist {output.allele_hist}
        --variants {output.variants}
        > {log} 2>&1
        """)


include: CFG["options"]["src"]["assets_workflow"]
include: CFG["options"]["src"]["regional_workflow"]


rule _aule_fragmentomics_output:
    input:
        summary = str(rules._aule_fragmentomics_collect.output.summary),
        size_hist = str(rules._aule_fragmentomics_collect.output.size_hist),
        motifs = str(rules._aule_fragmentomics_collect.output.motifs),
        allele_hist = str(rules._aule_fragmentomics_collect.output.allele_hist),
        variants = str(rules._aule_fragmentomics_collect.output.variants),
        regional_summary = str(rules._aule_fragmentomics_regional.output.summary),
        delfi = str(rules._aule_fragmentomics_regional.output.delfi),
        griffin_profile = str(rules._aule_fragmentomics_regional.output.griffin_profile),
        griffin_features = str(rules._aule_fragmentomics_regional.output.griffin_features)
    output:
        summary = CFG["dirs"]["outputs"] + "per_sample/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.fragmentomics.tsv",
        size_hist = CFG["dirs"]["outputs"] + "fragment_size/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.fragment_size_hist.tsv",
        motifs = CFG["dirs"]["outputs"] + "end_motifs/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.end_motifs.tsv",
        allele_hist = CFG["dirs"]["outputs"] + "allele_fragment_size/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.allele_fragment_size_hist.tsv",
        variants = CFG["dirs"]["outputs"] + "per_variant/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.per_variant.tsv",
        regional_summary = CFG["dirs"]["outputs"] + "per_sample/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.regional_fragmentomics.tsv",
        delfi = CFG["dirs"]["outputs"] + "delfi/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.delfi.tsv",
        griffin_profile = CFG["dirs"]["outputs"] + "griffin_profiles/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.griffin_profile.tsv",
        griffin_features = CFG["dirs"]["outputs"] + "griffin_features/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.griffin_features.tsv"
    run:
        for source, destination in zip(input, output):
            op.relative_symlink(source, destination, in_module=True)


rule _aule_fragmentomics_aggregate:
    input:
        summaries = expand(
            str(rules._aule_fragmentomics_collect.output.summary),
            zip,
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"]
        ),
        regional_summaries = expand(
            str(rules._aule_fragmentomics_regional.output.summary),
            zip,
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"]
        ),
        script = CFG["options"]["src"]["aggregate"]
    output:
        cohort = CFG["dirs"]["outputs"] + "aggregate/cohort_fragmentomics.tsv"
    params:
        inputs = lambda wildcards, input: ",".join(input.summaries),
        regional_inputs = lambda wildcards, input: ",".join(input.regional_summaries)
    conda:
        CFG["conda_envs"]["python"]
    resources:
        **CFG["resources"]["aggregate"]
    log:
        CFG["logs"]["outputs"] + "aggregate.log"
    shell:
        op.as_one_line("""
        python {input.script}
        --inputs '{params.inputs}' --regional-inputs '{params.regional_inputs}'
        --output {output.cohort}
        > {log} 2>&1
        """)


rule _aule_fragmentomics_all:
    input:
        expand(
            [
                str(rules._aule_fragmentomics_output.output.summary),
                str(rules._aule_fragmentomics_output.output.size_hist),
                str(rules._aule_fragmentomics_output.output.motifs),
                str(rules._aule_fragmentomics_output.output.allele_hist),
                str(rules._aule_fragmentomics_output.output.variants),
                str(rules._aule_fragmentomics_output.output.regional_summary),
                str(rules._aule_fragmentomics_output.output.delfi),
                str(rules._aule_fragmentomics_output.output.griffin_profile),
                str(rules._aule_fragmentomics_output.output.griffin_features)
            ],
            zip,
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"]
        ) + [str(rules._aule_fragmentomics_aggregate.output.cohort)]


op.cleanup_module(CFG)
