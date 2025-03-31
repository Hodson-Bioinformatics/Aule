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
# `CFG` is a shortcut to `config["lcr-modules"]["aule_umi_collapse"]`
CFG = op.setup_module(
    name = "aule_umi_collapse",
    version = "2.0",
    subdirectories = ["inputs", "keep_proper_pairs",
    "set_mate_info", "group_by_umi", "collapse_consensus_reads",
    "bwa", "merge_bam_alignment", "consensus_overlap", "sortbam", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _aule_umi_collapse_input_bam,
    _aule_umi_collapse_all

##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _aule_umi_collapse_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"],
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)

# Keep only properly paired reads - this step cannot be moved to the preproc stage
# because I suspect it may affect SVs calling when a mate is unmapped 
rule _aule_umi_collapse_keep_proper_pairs:
    input:
        bam = str(rules._aule_umi_collapse_input_bam.output.bam)
    output:
        bam = temp(CFG["dirs"]["keep_proper_pairs"] + "{seq_type}--{genome_build}/{sample_id}_properpairs.bam")
    log:
        stderr = CFG["logs"]["keep_proper_pairs"] + "{seq_type}--{genome_build}/{sample_id}.keep_proper_pairs.log"
    threads:
        CFG["threads"]["keep_proper_pairs"]
    conda:
        CFG["conda_envs"]["samtools"]
    shell:
        op.as_one_line("""
        samtools view -b -F 4,2048,256,512 -f 3 {input.bam} > {output.bam} 2> {log.stderr}
        """)
        
# Sort by query name 
rule _aule_umi_collapse_sortbam_queryname:
    input:
        bam = str(rules._aule_umi_collapse_keep_proper_pairs.output.bam)
    output:
        bam = temp(CFG["dirs"]["keep_proper_pairs"] + "{seq_type}--{genome_build}/{sample_id}_properpairs_byqueryname.bam")
    log:
        stdout = CFG["logs"]["keep_proper_pairs"] + "{seq_type}--{genome_build}/{sample_id}.sortbam_queryname.log",
        stderr = CFG["logs"]["keep_proper_pairs"] + "{seq_type}--{genome_build}/{sample_id}.sortbam_queryname.stderr.log"
    params:
        tempdir = CFG["scratch_directory"]
    resources: 
        **CFG["resources"]["sortbam_queryname"]
    threads:
        CFG["threads"]["sortbam_queryname"]
    conda:
        CFG["conda_envs"]["fgbio"]
    group: 
        "_aule_umi_collapse_fgbio_set_mate_{sample_id}"
    shell:
        op.as_one_line("""
        fgbio --tmp-dir={params.tempdir} {resources.mem_mb_java} 
        SortBam 
        -i {input.bam} 
        -o {output.bam} 
        -s Queryname
        > {log.stdout} 2> {log.stderr}
        """)

# Set read mate info 
rule _aule_umi_collapse_set_mate_info:
    input:
        bam = str(rules._aule_umi_collapse_sortbam_queryname.output.bam)
    output:
        bam = temp(CFG["dirs"]["set_mate_info"] + "{seq_type}--{genome_build}/{sample_id}_properpairs_mate.bam")
    log:
        stdout = CFG["logs"]["set_mate_info"] + "{seq_type}--{genome_build}/{sample_id}.set_mate_info.log",
        stderr = CFG["logs"]["set_mate_info"] + "{seq_type}--{genome_build}/{sample_id}.set_mate_info.stderr.log"
    params:
        tempdir = CFG["scratch_directory"]
    resources: 
        **CFG["resources"]["set_mate_info"]
    threads:
        CFG["threads"]["set_mate_info"]
    conda:
        CFG["conda_envs"]["fgbio"]
    group: 
        "_aule_umi_collapse_fgbio_set_mate_{sample_id}"
    shell:
        op.as_one_line("""
        fgbio --tmp-dir={params.tempdir} {resources.mem_mb_java}  
        SetMateInformation 
        -i {input.bam} 
        -o {output.bam} 
        --allow-missing-mates true
         > {log.stdout} 2> {log.stderr}
        """)

# Group by UMI
rule _aule_umi_collapse_group_by_umi:
    input:
        bam = str(rules._aule_umi_collapse_set_mate_info.output.bam)
    output:
        bam = temp(CFG["dirs"]["group_by_umi"] + "{seq_type}--{genome_build}/{sample_id}_groupedUMI.bam"),
        report = CFG["logs"]["group_by_umi"] + "{seq_type}--{genome_build}/{sample_id}_groupedUMI.txt"
    log:
        stdout = CFG["logs"]["group_by_umi"] + "{seq_type}--{genome_build}/{sample_id}.group_by_umi.log",
        stderr = CFG["logs"]["group_by_umi"] + "{seq_type}--{genome_build}/{sample_id}.group_by_umi.stderr.log"
    params:
        tempdir = CFG["scratch_directory"]
    resources: 
        **CFG["resources"]["group_by_umi"]
    threads:
        CFG["threads"]["group_by_umi"]
    conda:
        CFG["conda_envs"]["fgbio"]
    group: 
        "_aule_umi_collapse_fgbio_collapse_{sample_id}"
    shell:
        op.as_one_line("""
        fgbio --tmp-dir={params.tempdir} {resources.mem_mb_java} 
        GroupReadsByUmi 
        -i {input.bam} 
        -o {output.bam} 
        -t RX 
        --strategy adjacency 
        -f {output.report} 
        -l 9
        > {log.stdout} 2> {log.stderr}
        """)
 
# Sort by template coordinate 
rule _aule_umi_collapse_sort_coordinates:
    input:
        bam = str(rules._aule_umi_collapse_group_by_umi.output.bam)
    output:
        bam = temp(CFG["dirs"]["sortbam"] + "{seq_type}--{genome_build}/{sample_id}_groupedUMI_sorted.bam"),
    log:
        stdout = CFG["logs"]["sortbam"] + "{seq_type}--{genome_build}/{sample_id}.sort_coordinates.log",
        stderr = CFG["logs"]["sortbam"] + "{seq_type}--{genome_build}/{sample_id}.sort_coordinates.stderr.log"
    params:
        tempdir = CFG["scratch_directory"]
    resources: 
        **CFG["resources"]["sortbam_coordinates"]
    threads:
        CFG["threads"]["sortbam_coordinates"]
    conda:
        CFG["conda_envs"]["fgbio"]
    group: 
        "_aule_umi_collapse_fgbio_collapse_{sample_id}"
    shell:
        op.as_one_line("""
        fgbio --tmp-dir={params.tempdir} {resources.mem_mb_java} 
        SortBam 
        -i {input.bam} 
        -o {output.bam} 
        -s TemplateCoordinate
        > {log.stdout} 2> {log.stderr}
        """)
        
# Consensus reads collapse 
rule _aule_umi_collapse_consensus_reads:
    input:
        bam = str(rules._aule_umi_collapse_sort_coordinates.output.bam)
    output:
        ubam = temp(CFG["dirs"]["collapse_consensus_reads"] + "{seq_type}--{genome_build}/{sample_id}_consensus.ubam"),
    log:
        stdout = CFG["logs"]["collapse_consensus_reads"] + "{seq_type}--{genome_build}/{sample_id}.collapse_consensus_reads.log",
        stderr = CFG["logs"]["collapse_consensus_reads"] + "{seq_type}--{genome_build}/{sample_id}.collapse_consensus_reads.stderr.log"
    params:
        tempdir = CFG["scratch_directory"]
    resources: 
        **CFG["resources"]["consensus_reads"]
    threads:
        CFG["threads"]["consensus_reads"]
    conda:
        CFG["conda_envs"]["fgbio"]
    group: 
        "_aule_umi_collapse_fgbio_consensus_{sample_id}"
    shell:
        op.as_one_line("""
        fgbio --tmp-dir={params.tempdir} {resources.mem_mb_java} 
        CallMolecularConsensusReads 
        -i {input.bam} 
        -o {output.ubam} 
        --threads {threads} 
        --consensus-call-overlapping-bases false 
        --min-reads 1 -t MI 
        > {log.stdout} 2> {log.stderr}
        """)
        
# BWA input
rule _aule_umi_collapse_bwa_input:
    input:
        ubam = str(rules._aule_umi_collapse_consensus_reads.output.ubam)
    output:
        fastq = temp(CFG["dirs"]["bwa"] + "{seq_type}--{genome_build}/{sample_id}_bwa_input.fastq")
    log:
        stdout = CFG["logs"]["bwa"] + "{seq_type}--{genome_build}/{sample_id}.bwa_input.log",
        stderr = CFG["logs"]["bwa"] + "{seq_type}--{genome_build}/{sample_id}.bwa_input.stderr.log"
    params:
        options = CFG["options"]["bwa_input"]
    conda:
        CFG["conda_envs"]["picard"]
    resources: 
        **CFG["resources"]["bwa_mem"]
    threads:
        CFG["threads"]["bwa_mem"]
    group:
        "_aule_umi_collapse_bwa_{sample_id}"
    shell:
        op.as_one_line("""
        picard SamToFastq
        INPUT={input.ubam}
        FASTQ={output.fastq}
        {params.options}
         > {log.stdout} 2> {log.stderr}
        """)

# BWA MEM 
rule _aule_umi_collapse_bwa_mem:
    input:
        fastq = str(rules._aule_umi_collapse_bwa_input.output.fastq),
        ref = reference_files("genomes/{genome_build}/bwa_index/bwa-0.7.17/genome.fa")
    output:
        bam = temp(CFG["dirs"]["bwa"] + "{seq_type}--{genome_build}/{sample_id}_bwa.bam")
    log:
        stderr = CFG["logs"]["bwa"] + "{seq_type}--{genome_build}/{sample_id}_bwa.stderr.log"
    params:
        options = CFG["options"]["bwa_mem"]
    conda:
        CFG["conda_envs"]["bwa_aule"] 
    threads:
        CFG["threads"]["bwa_mem"]
    resources:
        **CFG["resources"]["bwa_mem"]
    group:
        "_aule_umi_collapse_bwa_{sample_id}"
    shell:
        op.as_one_line("""
        bwa mem {params.options}
        -t {threads} 
        {input.ref} 
        {input.fastq} | 
        samtools view -1 - > {output.bam}
        2> {log.stderr}
        """)
        
# Sort by query name (picard)
rule _aule_umi_collapse_ubam_sort_queryname:
    input:
        ubam = str(rules._aule_umi_collapse_consensus_reads.output.ubam)
    output:
        ubam = temp(CFG["dirs"]["bwa"] + "{seq_type}--{genome_build}/{sample_id}.sorted.ubam")
    log:
        stdout = CFG["logs"]["bwa"] + "{seq_type}--{genome_build}/{sample_id}.merge_bam_alignment.sort_queryname.log",
        stderr = CFG["logs"]["bwa"] + "{seq_type}--{genome_build}/{sample_id}.merge_bam_alignment.sort_queryname.stderr.log"
    resources: 
        **CFG["resources"]["bwa_mem"]
    conda:
        CFG["conda_envs"]["picard"]
    params:
        tempdir = CFG["scratch_directory"]
    group:
        "_aule_umi_collapse_bwa_{sample_id}"
    shell:
        op.as_one_line("""
        picard SortSam {resources.mem_mb_java} 
        --TMP_DIR {params.tempdir}
        -I {input.ubam} 
        -O {output.ubam} 
        -SO queryname 
        > {log.stdout} 2> {log.stderr}
        """)   
        
# Merge ubam and bam
rule _aule_umi_collapse_merge_bam_alignment:
    input:
        ubam = str(rules._aule_umi_collapse_ubam_sort_queryname.output.ubam),
        bam = str(rules._aule_umi_collapse_bwa_mem.output.bam),
        ref = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        bam = temp(CFG["dirs"]["merge_bam_alignment"] + "{seq_type}--{genome_build}/{sample_id}_collapsed.bam")
    log:
        stdout = CFG["logs"]["merge_bam_alignment"] + "{seq_type}--{genome_build}/{sample_id}.merge_bam_alignment.stdout.log",
        stderr = CFG["logs"]["merge_bam_alignment"] + "{seq_type}--{genome_build}/{sample_id}.merge_bam_alignment.stderr.log"
    params:
        options = CFG["options"]["merge_bam_alignment"],
        tempdir = CFG["scratch_directory"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["merge_bam_alignment"]
    resources:
        **CFG["resources"]["sortbam_queryname_picard"]
    group:
        "_aule_umi_collapse_postproc_{sample_id}"
    shell:
         op.as_one_line("""
        gatk MergeBamAlignment 
       --ALIGNED_BAM {input.bam} 
       --UNMAPPED_BAM {input.ubam} 
       --OUTPUT {output.bam} 
       --REFERENCE_SEQUENCE {input.ref} 
       --SORT_ORDER queryname
       {params.options} 
       --TMP_DIR {params.tempdir}
       > {log.stdout} 2> {log.stderr}
       """)
        
# Identify R1/R2 discordant bases
rule _aule_umi_collapse_consensus_overlap:
    input:
        bam = str(rules._aule_umi_collapse_merge_bam_alignment.output.bam),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        bam = temp(CFG["dirs"]["consensus_overlap"] + "{seq_type}--{genome_build}/{sample_id}.complete.bam"), 
        metrics = CFG["dirs"]["consensus_overlap"] + "{seq_type}--{genome_build}/{sample_id}.metrics.txt"
    log:
        stdout = CFG["logs"]["consensus_overlap"] + "{seq_type}--{genome_build}/{sample_id}.consensus_overlap.log",
        stderr = CFG["logs"]["consensus_overlap"] + "{seq_type}--{genome_build}/{sample_id}.consensus_overlap.log"
    params:
        tempdir = CFG["scratch_directory"]
    resources: 
        **CFG["resources"]["consensus_reads"]
    threads:
        CFG["threads"]["consensus_reads"]
    conda:
        CFG["conda_envs"]["fgbio"]
    shell:
        op.as_one_line("""
        samtools sort -n -u {input.bam} | 
        fgbio --tmp-dir={params.tempdir} {resources.mem_mb_java} 
        CallOverlappingConsensusBases 
        -i /dev/stdin
        -o {output.bam} 
        -m {output.metrics} 
        --threads {threads}
        -r {input.fasta}
        --sort-order Queryname
        --agreement-strategy PassThrough 
        --disagreement-strategy MaskBoth
        > {log.stdout} 2> {log.stderr}
        """)

# Sort by query name (picard)
rule _aule_umi_collapse_sort_final:
    input:
        bam = str(rules._aule_umi_collapse_consensus_overlap.output.bam)
    output:
        bam = CFG["dirs"]["sortbam"] + "{seq_type}--{genome_build}/{sample_id}.sorted.bam",
    log:
        stdout = CFG["logs"]["sortbam"] + "{seq_type}--{genome_build}/{sample_id}.merge_bam_alignment.sort_queryname.log",
        stderr = CFG["logs"]["sortbam"] + "{seq_type}--{genome_build}/{sample_id}.merge_bam_alignment.sort_queryname.stderr.log"
    params:
        tempdir = CFG["scratch_directory"]
    resources:
        **CFG["resources"]["sortbam_queryname_picard"]
    conda:
        CFG["conda_envs"]["picard"]
    group:
        "_aule_umi_collapse_postproc_{sample_id}"
    shell:
        op.as_one_line("""
        picard SortSam {resources.mem_mb_java}
        --TMP_DIR {params.tempdir}
        -I {input.bam}
        -O {output.bam}
        -SO queryname
        --CREATE_INDEX true
        > {log.stdout} 2> {log.stderr}
        """)
        
# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _aule_umi_collapse_output_bam:
    input:
        bam = str(rules._aule_umi_collapse_sort_final.output.bam),
    output:
        bam = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.bam",
    run:
        op.relative_symlink(input.bam, output.bam, in_module = True)

# Generates the target sentinels for each run, which generate the symlinks
rule _aule_umi_collapse_all:
    input:
        expand(
            [
                str(rules._aule_umi_collapse_output_bam.output.bam)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
