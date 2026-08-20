rule _aule_fragmentomics_delfi_bins:
    input:
        fasta = CFG["inputs"]["reference_fasta"],
        fai = CFG["inputs"]["reference_fasta_fai"],
        blacklist = CFG["inputs"]["blacklist_bed"],
        script = CFG["options"]["src"]["prepare_regional_assets"]
    output:
        DELFI_BINS_GENERATED
    params:
        bin_size = CFG["options"]["delfi_bin_size_bp"],
        chromosomes = CFG["options"]["delfi_chromosomes"],
        max_blacklist = CFG["options"]["delfi_max_blacklist_fraction"]
    conda:
        CFG["conda_envs"]["python"]
    threads:
        CFG["threads"]["assets"]
    resources:
        **CFG["resources"]["assets"]
    log:
        CFG["logs"]["outputs"] + "prepare_delfi_bins.log"
    shell:
        op.as_one_line("""
        python {input.script} delfi-bins
        --fasta {input.fasta}
        --fai {input.fai}
        --blacklist {input.blacklist}
        --bin-size {params.bin_size}
        --chromosomes '{params.chromosomes}'
        --max-blacklist-fraction {params.max_blacklist}
        --output {output}
        > {log} 2>&1
        """)


rule _aule_fragmentomics_griffin_sites:
    input:
        sites = _griffin_raw_site_dependencies,
        fai = CFG["inputs"]["reference_fasta_fai"],
        blacklist = CFG["inputs"]["blacklist_bed"],
        allowlist = _griffin_allowlist_dependencies,
        script = CFG["options"]["src"]["prepare_regional_assets"]
    output:
        GRIFFIN_SITES_GENERATED
    params:
        allowlist = lambda w: _formatted_input("griffin_tf_allowlist", w),
        window = CFG["options"]["griffin_window_bp"]
    conda:
        CFG["conda_envs"]["python"]
    threads:
        CFG["threads"]["assets"]
    resources:
        **CFG["resources"]["assets"]
    log:
        CFG["logs"]["outputs"] + "prepare_griffin_sites.log"
    shell:
        op.as_one_line("""
        python {input.script} griffin-sites
        --input {input.sites}
        --fai {input.fai}
        --blacklist {input.blacklist}
        --window {params.window}
        --tf-allowlist '{params.allowlist}'
        --output {output}
        > {log} 2>&1
        """)


rule _aule_fragmentomics_delfi_reference:
    input:
        controls = DELFI_REFERENCE_INPUTS,
        script = CFG["options"]["src"]["build_regional_reference"]
    output:
        DELFI_REFERENCE_GENERATED
    params:
        controls = lambda wildcards, input: ",".join(input.controls),
        min_samples = CFG["options"]["reference_min_samples"]
    conda:
        CFG["conda_envs"]["python"]
    threads:
        CFG["threads"]["reference"]
    resources:
        **CFG["resources"]["reference"]
    log:
        CFG["logs"]["outputs"] + "build_delfi_reference.log"
    shell:
        op.as_one_line("""
        python {input.script}
        --delfi-inputs '{params.controls}'
        --delfi-output {output}
        --min-samples {params.min_samples}
        > {log} 2>&1
        """)


rule _aule_fragmentomics_griffin_reference:
    input:
        controls = GRIFFIN_REFERENCE_INPUTS,
        script = CFG["options"]["src"]["build_regional_reference"]
    output:
        GRIFFIN_REFERENCE_GENERATED
    params:
        controls = lambda wildcards, input: ",".join(input.controls),
        min_samples = CFG["options"]["reference_min_samples"]
    conda:
        CFG["conda_envs"]["python"]
    threads:
        CFG["threads"]["reference"]
    resources:
        **CFG["resources"]["reference"]
    log:
        CFG["logs"]["outputs"] + "build_griffin_reference.log"
    shell:
        op.as_one_line("""
        python {input.script}
        --griffin-inputs '{params.controls}'
        --griffin-output {output}
        --min-samples {params.min_samples}
        > {log} 2>&1
        """)
