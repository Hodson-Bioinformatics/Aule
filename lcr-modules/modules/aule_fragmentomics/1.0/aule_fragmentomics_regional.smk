rule _aule_fragmentomics_regional:
    input:
        bam = str(rules._aule_fragmentomics_input_bam.output.bam),
        bai = str(rules._aule_fragmentomics_input_bam.output.bai),
        fasta = CFG["inputs"]["reference_fasta"],
        fasta_fai = CFG["inputs"]["reference_fasta_fai"],
        blacklist = CFG["inputs"]["blacklist_bed"],
        delfi_bins = lambda w: _selected_dependency(_delfi_bins, ENABLE_DELFI, w),
        delfi_reference = lambda w: _selected_dependency(_delfi_reference, ENABLE_DELFI, w),
        griffin_sites = lambda w: _selected_dependency(_griffin_sites, ENABLE_GRIFFIN, w),
        griffin_reference = lambda w: _selected_dependency(_griffin_reference, ENABLE_GRIFFIN, w),
        script = CFG["options"]["src"]["collect_regional"]
    output:
        summary = CFG["dirs"]["per_sample"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.regional_fragmentomics.tsv",
        delfi = CFG["dirs"]["per_sample"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.delfi.tsv",
        griffin_profile = CFG["dirs"]["per_sample"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.griffin_profile.tsv",
        griffin_features = CFG["dirs"]["per_sample"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.griffin_features.tsv"
    params:
        assay_partition = CFG["options"]["assay_partition"],
        enable_delfi = str(ENABLE_DELFI).lower(),
        enable_griffin = str(ENABLE_GRIFFIN).lower(),
        delfi_bins = lambda w: _delfi_bins(w) if ENABLE_DELFI else "",
        delfi_reference = lambda w: _delfi_reference(w) if ENABLE_DELFI else "",
        griffin_sites = lambda w: _griffin_sites(w) if ENABLE_GRIFFIN else "",
        griffin_reference = lambda w: _griffin_reference(w) if ENABLE_GRIFFIN else "",
        min_mapq = CFG["options"]["min_mapq"],
        exclude_flags = CFG["options"]["exclude_flags"],
        min_length = CFG["options"]["min_fragment_length"],
        max_length = CFG["options"]["max_fragment_length"],
        short_min = CFG["options"]["short_min"],
        short_max = CFG["options"]["short_max"],
        long_min = CFG["options"]["long_min"],
        long_max = CFG["options"]["long_max"],
        delfi_pseudocount = CFG["options"]["delfi_pseudocount"],
        delfi_min_fragments = CFG["options"]["delfi_min_fragments_per_bin"],
        delfi_min_bins = CFG["options"]["delfi_min_bins_for_lowess"],
        delfi_span = CFG["options"]["delfi_gc_lowess_span"],
        griffin_window = CFG["options"]["griffin_window_bp"],
        griffin_bin = CFG["options"]["griffin_bin_size_bp"],
        griffin_central = CFG["options"]["griffin_central_bp"],
        griffin_amplitude = CFG["options"]["griffin_amplitude_bp"],
        griffin_min_bins = CFG["options"]["griffin_min_bins_for_lowess"],
        griffin_min_sites = CFG["options"]["griffin_min_sites_per_tf"],
        griffin_span = CFG["options"]["griffin_gc_lowess_span"]
    conda:
        CFG["conda_envs"]["python"]
    threads:
        CFG["threads"]["regional"]
    resources:
        **CFG["resources"]["regional"]
    log:
        CFG["logs"]["per_sample"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.regional.log"
    shell:
        op.as_one_line("""
        python {input.script}
        --sample-id '{wildcards.tumour_id}'
        --assay-partition '{params.assay_partition}'
        --bam {input.bam}
        --fasta {input.fasta}
        --blacklist {input.blacklist}
        --enable-delfi {params.enable_delfi}
        --delfi-bins '{params.delfi_bins}'
        --delfi-reference '{params.delfi_reference}'
        --enable-griffin {params.enable_griffin}
        --griffin-sites '{params.griffin_sites}'
        --griffin-reference '{params.griffin_reference}'
        --min-mapq {params.min_mapq}
        --exclude-flags {params.exclude_flags}
        --min-length {params.min_length}
        --max-length {params.max_length}
        --short-min {params.short_min}
        --short-max {params.short_max}
        --long-min {params.long_min}
        --long-max {params.long_max}
        --delfi-pseudocount {params.delfi_pseudocount}
        --delfi-min-fragments {params.delfi_min_fragments}
        --delfi-min-bins-for-lowess {params.delfi_min_bins}
        --delfi-gc-span {params.delfi_span}
        --griffin-window-bp {params.griffin_window}
        --griffin-bin-size-bp {params.griffin_bin}
        --griffin-central-bp {params.griffin_central}
        --griffin-amplitude-bp {params.griffin_amplitude}
        --griffin-min-bins-for-lowess {params.griffin_min_bins}
        --griffin-min-sites {params.griffin_min_sites}
        --griffin-gc-span {params.griffin_span}
        --delfi-output {output.delfi}
        --griffin-profile-output {output.griffin_profile}
        --griffin-features-output {output.griffin_features}
        --summary-output {output.summary}
        > {log} 2>&1
        """)
