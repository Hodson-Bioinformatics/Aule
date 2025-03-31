#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Joanna A. Krupka
# Module Author:    Joanna A. Krupka
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op
import hashlib

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
# `CFG` is a shortcut to `config["lcr-modules"]["ichorcna_pon"]`
CFG = op.setup_module(
    name = "ichorcna_pon",
    version = "1.0",
    subdirectories = ["inputs", "bamCoverage", "bigwigToWig", "convert_wig", "compile_wigs", "ichorcna_pon", "outputs"],
)

# Define rules to be run locally when using a compute cluster=
localrules:
    _ichorcna_pon_input_bam,
    _ichorcna_pon_output_pon,
    _ichorcna_pon_all


##### RULES #####

# This defines the script/extdata directory used by ichorCNA in the subsequent rules:
# Obtain the path to the ichorCNA conda environment
md5hash = hashlib.md5()
if workflow.conda_prefix:
    conda_prefix = workflow.conda_prefix
else:
    conda_prefix = os.path.abspath(".snakemake/conda")

md5hash.update(conda_prefix.encode())
f = open(CFG["conda_envs"]["ichorcna_pon"], 'rb')
md5hash.update(f.read())
f.close()
h = md5hash.hexdigest()

ichorcna_env = conda_prefix + "/" + h + "_" 
ichorcna_dir = ichorcna_env + "/ichorCNA"
ichorcna_scripts_dir = conda_prefix + "/" + h + "_" + "/ichorCNA/scripts/"

# Install ichorCNA
rule _ichorcna_pon_install_ichorcna:
    params:
        branch = CFG['options']['ichorcna_branch'],
        directory = ichorcna_dir
    output:
        ichorcna = directory(ichorcna_env + "/lib/R/library/ichorCNA")
    log:
        stdout = CFG["logs"]["inputs"] + "install_ichorCNA.stdout.log",
        stderr = CFG["logs"]["inputs"] + "install_ichorCNA.stderr.log",
    conda:
        CFG['conda_envs']['ichorcna_pon']
    shell:
        op.as_one_line("""
        git clone --branch {params.branch} --single-branch https://github.com/broadinstitute/ichorCNA.git {params.directory} &&
        R -q -e "options(timeout=9999999); install.packages('{params.directory}', repos = NULL, type = 'source')"
        """)
        
# Symlinks the input files into the module results directory (under '00-inputs/')
rule _ichorcna_pon_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)

# Deeptools to get .bw from .bam
rule _ichorcna_pon_bamCoverage:
    input:
        bam = str(rules._ichorcna_pon_input_bam.output.bam)
    output:
        bw = temp(CFG["dirs"]["bamCoverage"] + "{seq_type}--{genome_build}/bin{binSize}kb/{sample_id}.bin{binSize}kb.bw")
    params:
        binSize = "{binSize}",
        opts = CFG["options"]["bamCoverage"]["opts"]
    conda: 
        CFG["conda_envs"]["bamCoverage"]
    threads: 
        CFG["threads"]["bamCoverage"]
    resources:
        **CFG["resources"]["bamCoverage"]
    wildcard_constraints: 
        binSize = "|".join(["10", "50", "500", "1000"])
    log:
        CFG["logs"]["bamCoverage"] + "{seq_type}--{genome_build}/bin{binSize}kb/{sample_id}.bin{binSize}kb.log"
    shell:
         """
             bamCoverage -b {input.bam} --binSize {params.binSize}000 {params.opts} -o {output.bw} -p {threads}
         """

# Converts bigWig to Wig
rule _ichorcna_pon_bigwigToWig:
    input:
        bw = str(rules._ichorcna_pon_bamCoverage.output.bw)
    output:
        wig = temp(CFG["dirs"]["bigwigToWig"] + "{seq_type}--{genome_build}/bin{binSize}kb/{sample_id}.bin{binSize}kb_{chrom}.wig")
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
        CFG["logs"]["bigwigToWig"] + "{seq_type}--{genome_build}/bin{binSize}kb/{sample_id}.bin{binSize}kb_{chrom}.log"
    shell:
        """
        bigWigToWig {input.bw} {output.wig} -chrom={wildcards.chrom} 
        """

# This function will reformat the wig file to one that can be used for ichorCNA
rule _ichorcna_pon_convert_wig:
    input:
        wig = str(rules._ichorcna_pon_bigwigToWig.output.wig)
    output:
        wig = temp(CFG["dirs"]["convert_wig"] + "{seq_type}--{genome_build}/bin{binSize}kb/{sample_id}.bin{binSize}kb.{chrom}--fixed.wig")
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
        CFG["logs"]["convert_wig"] + "{seq_type}--{genome_build}/bin{binSize}kb/{sample_id}.bin{binSize}kb.{chrom}--fixed.log"
    shell:
        op.as_one_line("""
             echo -e "fixedStep chrom={wildcards.chrom} start=1 step={wildcards.binSize}000 span={wildcards.binSize}000 " > {output.wig} &&
             intersectBed -a <( bedops --chop {wildcards.binSize}000 --header  {input.wig} ) -b {input.wig} -wa -wb | awk '{{print $7}}' >> {output.wig} 2>> {log}
         """)

# This function is used to get the wigs of the main chromosomes, which will be stitched together
def get_chrom_wigs(wildcards):
    CFG = config["lcr-modules"]["ichorcna_pon"]
    chrs = reference_files("genomes/" + wildcards.genome_build + "/genome_fasta/main_chromosomes_withY.txt")
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
        wigs = expand(
            CFG["dirs"]["convert_wig"] + "{{seq_type}}--{{genome_build}}/bin{{binSize}}kb/{{sample_id}}.bin{{binSize}}kb.{chrom}--fixed.wig", 
            chrom = chrs)
    return(wigs)

rule _ichorcna_pon_compile_wigs:
    input:
        wigs = get_chrom_wigs
    output:
        wig = CFG["dirs"]["compile_wigs"] + "{seq_type}--{genome_build}/bin{binSize}kb/{sample_id}.bin{binSize}kb.wig"
    shell:
      """
      cat {input.wigs} > {output.wig}
      """

# Generate a list of all files for the PON 
rule _ichorcna_pon_wigpaths:
    input:
        bw = expand(
            [
                str(rules._ichorcna_pon_compile_wigs.output.wig)
            ],
            zip,
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"],
            binSize= [CFG["options"]["binSize"]] * len(CFG["samples"]["sample_id"]),
            pon_suffix = CFG["options"]["pon_suffix"],
            allow_missing = True)
    output:
        txt = CFG["dirs"]["ichorcna_pon"] + "{seq_type}--{genome_build}/bin{binSize}kb/input_wig_list_{pon_suffix}.txt"
    shell:
        """
            readlink -f {input.bw} > {output.txt}
        """

# Generate PON 
def get_gcwig(wildcards):
    genome_ucsc = op.switch_on_wildcard("genome_build", config["lcr-modules"]["ichorcna_pon"]["options"]["ichorcna_pon"]["genome_ucsc"])(wildcards)
    return ichorcna_dir + "/inst/extdata/gc_" + genome_ucsc + "_" + wildcards.binSize + "kb.wig"

def get_mapwig(wildcards):
    genome_ucsc = op.switch_on_wildcard("genome_build", config["lcr-modules"]["ichorcna_pon"]["options"]["ichorcna_pon"]["genome_ucsc"])(wildcards)
    return ichorcna_dir + "/inst/extdata/map_" + genome_ucsc + "_" + wildcards.binSize + "kb.wig"

def get_centromere(wildcards):
    centromere = op.switch_on_wildcard("genome_build", config["lcr-modules"]["ichorcna_pon"]["options"]["ichorcna_pon"]["centromere"])(wildcards)
    return ichorcna_dir + "/inst/extdata/" + centromere

rule _ichorcna_pon_create_pon:
    input:
        wigs = str(rules._ichorcna_pon_wigpaths.output.txt),
        ichorCNA = str(rules._ichorcna_pon_install_ichorcna.output.ichorcna)
    output:
        pon_rds = CFG["dirs"]["ichorcna_pon"] + "{seq_type}--{genome_build}/bin{binSize}kb/ichorCNA_PON_{genome_build}_bin{binSize}_{pon_suffix}_{method}.rds",
        pon_txt = CFG["dirs"]["ichorcna_pon"] + "{seq_type}--{genome_build}/bin{binSize}kb/ichorCNA_PON_{genome_build}_bin{binSize}_{pon_suffix}_{method}.txt"
    params:
        ichorcna_scripts_dir = ichorcna_scripts_dir,
        genome_style = op.switch_on_wildcard("genome_build", CFG["options"]["ichorcna_pon"]["genome_style"]),
        gcwig = get_gcwig,
        mapwig = get_mapwig,
        centromere = get_centromere,
        prefix = CFG["dirs"]["ichorcna_pon"] + "{seq_type}--{genome_build}/bin{binSize}kb/ichorCNA_PON_{genome_build}_bin{binSize}_{pon_suffix}",
        opts = CFG["options"]["ichorcna_pon"]["opts"]
    conda: 
        CFG["conda_envs"]["ichorcna_pon"]
    threads: 
        CFG["threads"]["ichorcna_pon"]
    resources:
        **CFG["resources"]["ichorcna_pon"]
    log:
        CFG["logs"]["ichorcna_pon"] + "{seq_type}--{genome_build}/bin{binSize}kb/ichorcna_pon_{pon_suffix}_{method}.log"
    shell:
      op.as_one_line("""
      Rscript {params.ichorcna_scripts_dir}/createPanelOfNormals.R 
      --filelist {input.wigs} 
      -o {params.prefix} 
      --gcWig {params.gcwig} 
      --centromere {params.centromere}  
      --mapWig {params.mapwig} 
      {params.opts}
      """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _ichorcna_pon_output_pon:
    input:
        pon_rds = str(rules._ichorcna_pon_create_pon.output.pon_rds),
        pon_txt = str(rules._ichorcna_pon_create_pon.output.pon_txt)
    output:
        pon_rds = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/bin{binSize}kb/ichorCNA_PON_{genome_build}_bin{binSize}_{pon_suffix}_{method}.rds",
        pon_txt = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/bin{binSize}kb/ichorCNA_PON_{genome_build}_bin{binSize}_{pon_suffix}_{method}.txt"
    run:
        op.relative_symlink(input.pon_rds, output.pon_rds, in_module= True)
        op.relative_symlink(input.pon_txt, output.pon_txt, in_module= True)

# Generates the target sentinels for each run, which generate the symlinks
rule _ichorcna_pon_all:
    input:
        expand(
            [
                str(rules._ichorcna_pon_output_pon.output.pon_rds),
                str(rules._ichorcna_pon_output_pon.output.pon_txt)
            ],
            zip,
            seq_type = CFG["samples"]["seq_type"],
            genome_build = CFG["samples"]["genome_build"],
            binSize = CFG["options"]["binSize"],
            method = CFG["options"]["ichorcna_pon"]["method"],
            pon_suffix = CFG["options"]["pon_suffix"],
            allow_missing = True
        )


##### CLEANUP #####

# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
