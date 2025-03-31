#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Joanna A. Krupka 
# Module Author:    Joanna A. Krupka
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op
import pandas as pd

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
# `CFG` is a shortcut to `config["lcr-modules"]["aule_pv"]`
CFG = op.setup_module(
    name = "aule_pv",
    version = "1.0",
    subdirectories = ["inputs", "ref", "pv_all", "sorted_bam", "pv_collapsed", "annotated", "genotype", "track", "mrd", "report","outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _aule_pv_input_bam,
    _aule_pv_input_bed,
    _aule_pv_output_txt,
    _aule_pv_all

# Define PV output
pv_runs = pd.concat([CFG["runs"], CFG["runs"]])
pv_runs["pv_type"] = ["doublets"] * len(CFG["runs"]) + ["triplets"] * len(CFG["runs"])
pv_runs['tumour_patient_id'] = pv_runs['tumour_patient_id'].str.split('_').str[0]
pv_out = pv_runs.loc[:, ["tumour_patient_id", "tumour_seq_type", "tumour_genome_build", "pv_type"]].drop_duplicates()

##### RULES #####
# Symlinks the input files into the module results directory (under '00-inputs/')
rule _aule_pv_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        op.absolute_symlink(input.bam, output.bam)

rule _aule_pv_input_bed:
    input:
        bed = CFG["inputs"]["panel"]
    output:
        bed = CFG["dirs"]["inputs"] + "bed/target_regions.bed"
    run:
        op.absolute_symlink(input.bed, output.bed)
        
# Temporary reference rule (awaiting permission change)
rule chromsizes_genome_fasta:
    input: 
        fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai")
    output: 
        chromsizes = CFG["dirs"]["ref"] + "genomes/{genome_build}/genome.fa.chromsizes"
    conda:
        CFG["conda_envs"]["coreutils"]
    shell:
        "cut -f1,2 {input.fai} > {output.chromsizes}"
        
# Find phased variants in tumour samples 
rule _aule_pv_find_all_phased_variants:
    input:
        bam = str(rules._aule_pv_input_bam.output.bam) ,
        chromsizes = str(rules.chromsizes_genome_fasta.output.chromsizes),
    output:
        doublets = CFG["dirs"]["pv_all"] + "{seq_type}--{genome_build}/{sample_id}.doublets.all.txt",
        triplets = CFG["dirs"]["pv_all"] + "{seq_type}--{genome_build}/{sample_id}.triplets.all.txt",
        stat = CFG["dirs"]["pv_all"] + "{seq_type}--{genome_build}/{sample_id}.pv_all_stats.txt",
    log:
        stdout = CFG["logs"]["pv_all"] + "{seq_type}--{genome_build}/{sample_id}.doublets.find_all_pv.stdout.log",
        stderr = CFG["logs"]["pv_all"] + "{seq_type}--{genome_build}/{sample_id}.triplets.find_all_pv.stderr.log"
    params:
        script = CFG["scripts"]["find_all_phased_variants"],
        min_mapping_qual = CFG["options"]["min_mapping_qual"],
        min_base_qual = CFG["options"]["min_base_qual"]
    conda:
        CFG["conda_envs"]["phased_variants"]
    threads:
        CFG["threads"]["find_all_phased_variants"]
    resources:
        **CFG["resources"]["find_all_phased_variants"]
    shell:
        op.as_one_line("""
        python {params.script}
               -bam {input.bam}
               -chromsizes {input.chromsizes}
               -output_doublets {output.doublets}
               -output_triplets {output.triplets}
               -stat {output.stat}
               -min_mapping_qual {params.min_mapping_qual}
               -min_base_qual {params.min_base_qual}
               > {log.stdout} 2> {log.stderr}
        """)

# Sort by coordinates
rule _aule_pv_sort_by_coordinate: 
    input:  
        bam = str(rules._aule_pv_input_bam.output.bam),
    output:
        bam = CFG["dirs"]["sorted_bam"] + "{seq_type}--{genome_build}/{sample_id}.sorted_coordinate.bam"
    log: 
        stdout = CFG["logs"]["sorted_bam"] + "{seq_type}--{genome_build}/{sample_id}.sorted_coordinate.stdout.log",
        stderr = CFG["logs"]["sorted_bam"] + "{seq_type}--{genome_build}/{sample_id}.sorted_coordinate.stderr.log"
    params:
        index =  CFG["dirs"]["sorted_bam"] + "{seq_type}--{genome_build}/{sample_id}.sorted_coordinate.bai",
        tempdir = CFG["scratch_directory"]
    conda: 
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["sort_bam"]
    resources:
        **CFG["resources"]["sort_bam"]
    shell:
        op.as_one_line("""
        gatk SortSam 
        --INPUT {input.bam} 
        --OUTPUT {output.bam} 
        --SORT_ORDER 'coordinate' 
        --CREATE_INDEX true 
        --CREATE_MD5_FILE false
        --TMP_DIR {params.tempdir}
        && touch {params.index}.bai
         """)

# Collapse and filter
rule _aule_pv_collapse_and_filter: 
    input:
        doublets = str(rules._aule_pv_find_all_phased_variants.output.doublets),
        triplets = str(rules._aule_pv_find_all_phased_variants.output.triplets)
    output:
        doublets = CFG["dirs"]["pv_collapsed"] + "{seq_type}--{genome_build}/{sample_id}.doublets.raw.txt",    
        triplets = CFG["dirs"]["pv_collapsed"] + "{seq_type}--{genome_build}/{sample_id}.triplets.raw.txt",
        bed = CFG["dirs"]["pv_collapsed"] + "{seq_type}--{genome_build}/{sample_id}.pv_all.bed"
    log:
        stdout = CFG["logs"]["pv_collapsed"] + "{seq_type}--{genome_build}/{sample_id}.pv_collapsed.stdout.log",
        stderr = CFG["logs"]["pv_collapsed"] + "{seq_type}--{genome_build}/{sample_id}.pv_collapsed.stderr.log"
    params:
        min_mapping_qual = CFG["options"]["min_mapping_qual"],
        min_base_qual = CFG["options"]["min_base_qual"],
        r1r2 = CFG["options"]["r1r2"],
        detailed = CFG["options"]["detailed"]
    conda:
        CFG["conda_envs"]["postprocessing"]
    threads:
        CFG["threads"]["collapse_phased_variants"]
    resources:
        **CFG["resources"]["collapse_phased_variants"]
    script:
        "etc/R/PV_collapse_and_filter.R"

# Annotate PVs with matched germline, 3nt context, total DP and and read support in healthy cohort (PON)
rule _aule_pv_annotate:
    input:
        tumour_bam = CFG["dirs"]["sorted_bam"] + "{seq_type}--{genome_build}/{tumour_id}.sorted_coordinate.bam",
        normal_bam = CFG["dirs"]["sorted_bam"] + "{seq_type}--{genome_build}/{normal_id}.sorted_coordinate.bam",
        tumour_pv = CFG["dirs"]["pv_collapsed"] + "{seq_type}--{genome_build}/{tumour_id}.{pv_type}.raw.txt",
        normal_pv = CFG["dirs"]["pv_collapsed"] + "{seq_type}--{genome_build}/{normal_id}.{pv_type}.raw.txt",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        panel = str(rules._aule_pv_input_bed.output.bed),
        pon_pv = CFG["inputs"]["pon_pv"],
        pon_snv = CFG["inputs"]["pon_snv"], 
        gnomad = CFG["inputs"]["gnomad"]
    output:
        pv_annotated = CFG["dirs"]["annotated"] + "{seq_type}--{genome_build}/{tumour_id}-{normal_id}.{pv_type}.annotated.txt",
    log:
        stdout = CFG["logs"]["annotated"] + "{seq_type}--{genome_build}/{tumour_id}-{normal_id}.{pv_type}.annotated.stdout.log",
        stderr = CFG["logs"]["annotated"] + "{seq_type}--{genome_build}/{tumour_id}-{normal_id}.{pv_type}.annotated.stderr.log",
    params:
        script = CFG["scripts"]["annotate_phased_variants"],
        pv_type = "{pv_type}"
    conda:
        CFG["conda_envs"]["phased_variants"]
    threads:
        CFG["threads"]["annotate_phased_variants"]
    resources:
        **CFG["resources"]["annotate_phased_variants"]
    shell:
        op.as_one_line("""
        python {params.script} 
        -tumour_bam {input.tumour_bam}
        -normal_bam {input.normal_bam}
        -tumour_pv {input.tumour_pv} 
        -normal_pv {input.normal_pv}
        -pv_type {params.pv_type}
        -fasta {input.fasta}
        -pon_pv {input.pon_pv}
        -pon_snv {input.pon_snv}
        -gnomad {input.gnomad}
        -panel {input.panel}
        -pv_out {output.pv_annotated} 
        """)

# Perform genotyping 
rule _aule_pv_genotype:
    input:
        pv_annotated = str(rules._aule_pv_annotate.output.pv_annotated),
    output:
        pv_genotyped = CFG["dirs"]["genotype"] + "{seq_type}--{genome_build}/{tumour_id}-{normal_id}.{pv_type}.genotype.txt",
    log:
        stdout = CFG["logs"]["genotype"] + "{seq_type}--{genome_build}/{tumour_id}-{normal_id}.{pv_type}.genotype.stdout.log",
        stderr = CFG["logs"]["genotype"] + "{seq_type}--{genome_build}/{tumour_id}-{normal_id}.{pv_type}.genotype.stderr.log",
    params:
        pv_type = "{pv_type}",
        max_matched_germline_support = CFG["options"]["genotype_PV"]["max_matched_germline_support"],
        max_cohort_germline_support = CFG["options"]["genotype_PV"]["max_cohort_germline_support"],
        max_gnomad_support = CFG["options"]["genotype_PV"]["max_gnomad_support"],
        max_pon_snv_support = CFG["options"]["genotype_PV"]["max_pon_snv_support"],
        max_pon_pv_support = CFG["options"]["genotype_PV"]["max_pon_pv_support"],
        min_AF = CFG["options"]["genotype_PV"]["min_AF"],
        min_DP = CFG["options"]["genotype_PV"]["min_DP"],
    conda:
        CFG["conda_envs"]["postprocessing"]
    threads:
        CFG["threads"]["genotype_phased_variants"]
    resources:
        **CFG["resources"]["genotype_phased_variants"]
    script:
        "etc/R/PV_genotype.R"
    
# Track phased variants 
def get_pv_by_patient(wildcards):
    CFG = config["lcr-modules"]["aule_pv"]
    # Retreive all samples matching this patient
    this_patient = op.filter_samples(pv_runs, tumour_patient_id = wildcards.patient_id, tumour_genome_build = wildcards.genome_build, pv_type = wildcards.pv_type)
    
    pv_genotypes = expand(
            [
                CFG["dirs"]["genotype"] + "{seq_type}--{genome_build}/{tumour_id}-{normal_id}.{pv_type}.genotype.txt"
            ],
            zip,
            seq_type=this_patient["tumour_seq_type"],
            genome_build=this_patient["tumour_genome_build"],
            tumour_id=this_patient["tumour_sample_id"],
            normal_id=this_patient["normal_sample_id"],
            pair_status=this_patient["pair_status"],
            pv_type=this_patient["pv_type"],
            allow_missing = True)
            
    pv_all = expand(
            [
                CFG["dirs"]["pv_all"] + "{seq_type}--{genome_build}/{tumour_id}.{pv_type}.all.txt"
            ],
            zip,
            seq_type=this_patient["tumour_seq_type"],
            genome_build=this_patient["tumour_genome_build"],
            tumour_id=this_patient["tumour_sample_id"],
            pair_status=this_patient["pair_status"],
            pv_type=this_patient["pv_type"],
            allow_missing = True)
            
    bams = expand(
            [
                CFG["dirs"]["sorted_bam"] + "{seq_type}--{genome_build}/{tumour_id}.sorted_coordinate.bam"
            ],
            zip,
            seq_type=this_patient["tumour_seq_type"],
            genome_build=this_patient["tumour_genome_build"],
            tumour_id=this_patient["tumour_sample_id"],
            pair_status=this_patient["pair_status"],
            pv_type=this_patient["pv_type"],
            allow_missing = True)

    return {"pv_genotypes": pv_genotypes, "pv_all": pv_all, "bams": bams}

# Identify PVs to track (from baseline biopsy, baseline plasma or progressed plasma/biopsy) 
# Compute informative depth, PVAF per sample
rule _aule_pv_track:
    input:
        unpack(get_pv_by_patient)
    output: 
        pv_tracked = CFG["dirs"]["track"] + "{seq_type}--{genome_build}/{patient_id}.{pv_type}.tracked_pv.RData",
        pv_summarised = CFG["dirs"]["track"] + "{seq_type}--{genome_build}/{patient_id}.{pv_type}.summarised_pv.txt",
    log:
        stdout = CFG["logs"]["track"] + "{seq_type}--{genome_build}/{patient_id}.{pv_type}.tracked.stdout.log",
        stderr = CFG["logs"]["track"] + "{seq_type}--{genome_build}/{patient_id}.{pv_type}.tracked.stderr.log",
    params:
        pv_type = "{pv_type}",
        min_mapping_qual = CFG["options"]["min_mapping_qual"],
        min_base_qual = CFG["options"]["min_base_qual"],
    conda:
        CFG["conda_envs"]["postprocessing"]
    threads:
        CFG["threads"]["track_phased_variants"]
    resources:
        **CFG["resources"]["track_phased_variants"]
    script:
        "etc/R/PV_track.R"

# Compute MRD statistics: Monte Carlo simulation for signifcance of phased variants presence 
rule _aule_pv_mrd_statistic:
    input:
        unpack(get_pv_by_patient),
        pv_tracked = str(rules._aule_pv_track.output.pv_tracked),
        panel = str(rules._aule_pv_input_bed.output.bed)
    output: 
        pv_mrd = CFG["dirs"]["mrd"] + "{seq_type}--{genome_build}/{patient_id}.{pv_type}.mrd.txt",
        pv_tab = CFG["dirs"]["mrd"] + "{seq_type}--{genome_build}/{patient_id}.{pv_type}.mrd.pv.RData",
    log:
        stdout = CFG["logs"]["mrd"] + "{seq_type}--{genome_build}/{patient_id}.{pv_type}.mrd.stdout.log",
        stderr = CFG["logs"]["mrd"] + "{seq_type}--{genome_build}/{patient_id}.{pv_type}.mrd.stderr.log",
    params:
        pv_type = "{pv_type}",
        monte_carlo_iterations = CFG["options"]["mrd_stat"]["monte_carlo_iterations"],
        min_mapping_qual = CFG["options"]["min_mapping_qual"]
    conda:
        CFG["conda_envs"]["postprocessing"]
    threads:
        CFG["threads"]["mrd_phased_variants"]
    resources:
        **CFG["resources"]["mrd_phased_variants"]
    script:
        "etc/R/PV_stat.R"
        
    
# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _aule_pv_output_txt:
    input:
        pv_mrd = str(rules._aule_pv_mrd_statistic.output.pv_mrd)
    output:
        pv_mrd = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{patient_id}.{pv_type}.mrd.txt",
    run:
        op.relative_symlink(input.pv_mrd, output.pv_mrd, in_module= True)

# Output
rule _aule_pv_all:
    input:
        expand(
            [
                str(rules._aule_pv_output_txt.output.pv_mrd),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=pv_out["tumour_seq_type"],
            genome_build=pv_out["tumour_genome_build"],
            patient_id=pv_out["tumour_patient_id"],
            pv_type=pv_out["pv_type"],
            allow_missing=True
            )
            
##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
