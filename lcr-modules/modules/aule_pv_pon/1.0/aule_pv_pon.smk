#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Joanna A. Krupka
# Module Author:    Joanna A. Krupka
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
# `CFG` is a shortcut to `config["lcr-modules"]["aule_pv_pon"]`
CFG = op.setup_module(
    name = "aule_pv_pon",
    version = "1.0",
    subdirectories = ["inputs", "sorted_bam", "pv_all", "pv_collapsed", "pv_annotated","pv_pon", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _aule_pv_pon_input_bam,
    _aule_pv_pon_input_bed,
    _aule_pv_pon_chromsizes_genome_fasta,
    _aule_pv_all

##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _aule_pv_pon_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        op.absolute_symlink(input.bam, output.bam)

rule _aule_pv_pon_input_bed:
    input:
        bed = CFG["inputs"]["panel"]
    output:
        bed = CFG["dirs"]["inputs"] + "bed/target_regions.bed"
    run:
        op.absolute_symlink(input.bed, output.bed)
        
# Temporary reference rule
rule _aule_pv_pon_chromsizes_genome_fasta:
    input: 
        fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai")
    output: 
        chromsizes = CFG["dirs"]["inputs"] + "genomes/{genome_build}/genome.fa.chromsizes"
    conda:
        CFG["conda_envs"]["coreutils"]
    shell:
        "cut -f1,2 {input.fai} > {output.chromsizes}"

# Find phased variants in tumour samples 
rule _aule_pv_pon_find_all_phased_variants:
    input:
        bam = str(rules._aule_pv_pon_input_bam.output.bam),
        chromsizes = str(rules._aule_pv_pon_chromsizes_genome_fasta.output.chromsizes),
    output:
        doublets = temp(CFG["dirs"]["pv_all"] + "{seq_type}--{genome_build}/{sample_id}.doublets.all.txt"),
        triplets = temp(CFG["dirs"]["pv_all"] + "{seq_type}--{genome_build}/{sample_id}.triplets.all.txt"),
        bed = temp(CFG["dirs"]["pv_all"] + "{seq_type}--{genome_build}/{sample_id}.pv_all_regions.bed"),
        stat = temp(CFG["dirs"]["pv_all"] + "{seq_type}--{genome_build}/{sample_id}.pv_all_stats.txt"),
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
               -bed {output.bed}
               -min_mapping_qual {params.min_mapping_qual}
               -min_base_qual {params.min_base_qual}
               > {log.stdout} 2> {log.stderr}
        """)
    
# Sort by coordinates
rule _aule_pv_pon_sort_by_coordinate: 
    input:  
        bam = str(rules._aule_pv_pon_input_bam.output.bam),
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
rule _aule_pv_pon_collapse_and_filter: 
    input:
        doublets = str(rules._aule_pv_pon_find_all_phased_variants.output.doublets),
        triplets = str(rules._aule_pv_pon_find_all_phased_variants.output.triplets),
        bam = str(rules._aule_pv_pon_input_bam.output.bam),
    output:
        doublets = CFG["dirs"]["pv_collapsed"] + "{seq_type}--{genome_build}/{sample_id}.doublets.raw.txt",    
        triplets = CFG["dirs"]["pv_collapsed"] + "{seq_type}--{genome_build}/{sample_id}.triplets.raw.txt",
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
        
rule _aule_pv_pon_annotate: 
    input:
        doublets = str(rules._aule_pv_pon_collapse_and_filter.output.doublets),
        triplets = str(rules._aule_pv_pon_collapse_and_filter.output.triplets),
        bam = str(rules._aule_pv_pon_sort_by_coordinate.output.bam),
        panel = str(rules._aule_pv_pon_input_bed.output.bed)
    output:
        doublets = CFG["dirs"]["pv_annotated"] + "{seq_type}--{genome_build}/{sample_id}.doublets.annotated.txt",    
        triplets = CFG["dirs"]["pv_annotated"] + "{seq_type}--{genome_build}/{sample_id}.triplets.annotated.txt",
    log:
        stdout = CFG["logs"]["pv_annotated"] + "{seq_type}--{genome_build}/{sample_id}.annotated.log",
        stderr = CFG["logs"]["pv_annotated"] + "{seq_type}--{genome_build}/{sample_id}.annotated.log"
    params:
        script = CFG["scripts"]["annotate_phased_variants"],
        min_mapping_qual = CFG["options"]["min_mapping_qual"],
        min_base_qual = CFG["options"]["min_base_qual"]
    conda:
        CFG["conda_envs"]["phased_variants"]
    threads:
        CFG["threads"]["annotate_phased_variants"]
    resources:
        **CFG["resources"]["annotate_phased_variants"]
    shell:
        op.as_one_line("""
        python {params.script}
               -bam {input.bam}
               -input_doublets {input.doublets}
               -input_triplets {input.triplets}
               -output_doublets {output.doublets}
               -output_triplets {output.triplets}
               -panel {input.panel}
               -min_mapping_qual {params.min_mapping_qual}
               -min_base_qual {params.min_base_qual}
               > {log.stdout} 2> {log.stderr}
        """)
  
def get_pvs(wildcards):
    CFG = config["lcr-modules"]["aule_pv_pon"]
    all_samples = CFG["samples"]
    
    doublets = expand(
            [
                CFG["dirs"]["pv_annotated"] + "{seq_type}--{genome_build}/{sample_id}.doublets.annotated.txt"
            ],
            zip,
            seq_type=all_samples["seq_type"],
            genome_build=all_samples["genome_build"],
            sample_id=all_samples["sample_id"],
            allow_missing = True)
            
    triplets = expand(
            [
                CFG["dirs"]["pv_annotated"] + "{seq_type}--{genome_build}/{sample_id}.triplets.annotated.txt"
            ],
            zip,
            seq_type=all_samples["seq_type"],
            genome_build=all_samples["genome_build"],
            sample_id=all_samples["sample_id"],
            allow_missing = True)
            
    return({"doublets": doublets, "triplets": triplets})
 
    
# Generate combined PON
rule _aule_pv_pon_generate_pon: 
    input:
        unpack(get_pvs),
    output:
        pon_doublets = CFG["dirs"]["pv_pon"] + "{seq_type}--{genome_build}/pon_doublets.txt", 
        pon_triplets = CFG["dirs"]["pv_pon"] + "{seq_type}--{genome_build}/pon_triplets.txt"
    params:
        pon_samples = CFG["options"]["pon_samples"]
    conda:
        CFG["conda_envs"]["postprocessing"]
    threads:
        CFG["threads"]["generate_pon"]
    resources:
        **CFG["resources"]["generate_pon"]
    script:
        "etc/R/PV_pon.R"

# Output function

# Generates the target sentinels for each run, which generate the symlinks
pon_out = CFG["samples"].loc[:, ["seq_type", "genome_build"]].drop_duplicates()

rule _aule_pv_pon_all:
    input:
        expand(
            [
                str(rules._aule_pv_pon_generate_pon.output.pon_doublets),
                str(rules._aule_pv_pon_generate_pon.output.pon_triplets)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=pon_out["seq_type"],
            genome_build=pon_out["genome_build"], 
            allow_missing = True)


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
