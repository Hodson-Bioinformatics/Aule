#!/usr/bin/env snakemake

# ichorcna_burden_qc: ploidy-free copy-number burden index and noise metrics,
# computed from an EXISTING ichorcna_offtarget fit.
#
# This module never computes coverage and never runs ichorCNA. It reads the
# corrected depth and segmentation that ichorcna_offtarget published under its
# 99-outputs/, so every burden metric describes that same fit. Set
# inputs.ichorcna_output_dir to that 99-outputs root; the module stops at
# startup if it is unset.

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
    name = "ichorcna_burden_qc",
    version = "1.0",
    subdirectories = [
        "inputs",
        "burden",
        "burden_plot",
        "outputs"
    ]
)

localrules:
    _ichorcna_burden_qc_input_pon,
    _ichorcna_burden_qc_output,
    _ichorcna_burden_qc_aggregate,
    _ichorcna_burden_qc_all


##### ichorcna_offtarget OUTPUT #####

_ICHORCNA_OUTPUT_DIR = CFG["inputs"].get("ichorcna_output_dir", "")
if not _ICHORCNA_OUTPUT_DIR or _ICHORCNA_OUTPUT_DIR == "__UPDATE__":
    sys.exit(
        "\nichorcna_burden_qc: set inputs.ichorcna_output_dir to the 99-outputs\n"
        "directory of an ichorcna_offtarget run, e.g.\n"
        "    <results>/ichorcna_offtarget-1.0/99-outputs\n"
        "This module does not run ichorCNA itself.\n"
    )

# Layout of ichorcna_offtarget's 99-outputs/: subdirectory and file suffix.
_ICHOR_OUTPUT_LAYOUT = {
    "corrDepth": ("corrDepth", ".corrDepth.txt"),
    "segTxt": ("seg_txt", ".seg.txt"),
}


def _get_ichorcna_output(wildcards, output_name):
    output_directory, suffix = _ICHOR_OUTPUT_LAYOUT[output_name]
    return os.path.join(
        _ICHORCNA_OUTPUT_DIR.format(**wildcards),
        f"{wildcards.seq_type}--{wildcards.genome_build}",
        output_directory,
        f"bin{wildcards.binSize}kb",
        f"{wildcards.tumour_id}{suffix}"
    )


_get_burden_corr_depth = lambda wildcards: _get_ichorcna_output(
    wildcards, "corrDepth"
)
_get_burden_seg_txt = lambda wildcards: _get_ichorcna_output(
    wildcards, "segTxt"
)


##### RULES #####

# The PON is not used to normalise anything here (ichorcna_offtarget already
# did that); it is linked in so the burden table can record which PON the fit
# was made with. It must be the same PON the fit used.
rule _ichorcna_burden_qc_input_pon:
    input:
        pon = CFG["inputs"]["ichorcna_pon"]
    output:
        pon = CFG["dirs"]["inputs"] + "pon/{seq_type}--{genome_build}/bin{binSize}kb/ichorCNA_PON_{genome_build}_bin{binSize}_full_median.rds"
    log:
        CFG["logs"]["inputs"] + "pon/{seq_type}--{genome_build}/bin{binSize}kb.log"
    run:
        op.absolute_symlink(input.pon, output.pon)


rule _ichorcna_burden_qc_burden:
    input:
        bins = _get_burden_corr_depth,
        seg = _get_burden_seg_txt,
        pon = str(rules._ichorcna_burden_qc_input_pon.output.pon),
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
        CFG["conda_envs"]["r"]
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


rule _ichorcna_burden_qc_plot_vc_diagnostic:
    input:
        table = str(rules._ichorcna_burden_qc_burden.output.tsv),
        script = CFG["options"]["src"]["plot_vc_diagnostic"]
    output:
        pdf = CFG["dirs"]["burden_plot"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.vc_diagnostic.pdf"
    params:
        sample = "{tumour_id}"
    conda:
        CFG["conda_envs"]["r"]
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


# Symlinks the final output files into the module results directory (under '99-outputs/').
# The ichorCNA files themselves stay in ichorcna_offtarget's 99-outputs/.
rule _ichorcna_burden_qc_output:
    input:
        burden = str(rules._ichorcna_burden_qc_burden.output.tsv),
        burden_diagnostic = str(rules._ichorcna_burden_qc_plot_vc_diagnostic.output.pdf)
    output:
        burden = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/burden/bin{binSize}kb/{tumour_id}.burden.tsv",
        burden_diagnostic = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/burden_diagnostic/bin{binSize}kb/{tumour_id}.vc_diagnostic.pdf"
    log:
        CFG["logs"]["outputs"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.log"
    run:
        op.relative_symlink(input.burden, output.burden, in_module=True)
        op.relative_symlink(input.burden_diagnostic, output.burden_diagnostic, in_module=True)


rule _ichorcna_burden_qc_aggregate:
    input:
        burdens = expand(
            str(rules._ichorcna_burden_qc_output.output.burden),
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
        CFG["conda_envs"]["r"]
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


rule _ichorcna_burden_qc_all:
    input:
        expand(
            [
                str(rules._ichorcna_burden_qc_output.output.burden),
                str(rules._ichorcna_burden_qc_output.output.burden_diagnostic)
            ],
            zip,
            seq_type = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id = CFG["runs"]["tumour_sample_id"],
            binSize = [CFG["options"]["binsize_kb"]] * len(CFG["runs"]["tumour_sample_id"])
        ) + [
            str(rules._ichorcna_burden_qc_aggregate.output.cohort),
            str(rules._ichorcna_burden_qc_aggregate.output.qc)
        ]


##### CLEANUP #####

# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
