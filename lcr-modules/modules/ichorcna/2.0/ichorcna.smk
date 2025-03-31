#!/usr/bin/env snakemake


# ---------------------------------------------------------------------------- #
##### ATTRIBUTION #####
# ---------------------------------------------------------------------------- #

# Original snakemake author: Jasper Wong
# Module author: Ilias Moutsopoulos & Joanna Krupka & Rebecca Goodhew
# Additional contributors: N/A


# ---------------------------------------------------------------------------- #
##### SETUP #####
# ---------------------------------------------------------------------------- #

### Modules ###

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
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section 

### Directories ###
# Setup module and store module-specific configuration in `CFG`.
CFG = op.setup_module(
    name = "ichorcna", 
    version = "2.0",
    subdirectories = ["inputs", "bamCoverage", "bigwigToWig", "convert_wig", "compile_wigs", "ichorcna_run", "outputs"]
)

localrules:
    _ichorcna_input_bam,
    _ichorcna_input_pon,
    _ichorcna_output,
    _ichorcna_all

# ---------------------------------------------------------------------------- #
##### RULES #####
# ---------------------------------------------------------------------------- #

# This defines the script/extdata directory used by ichorCNA in the subsequent rules:
# Obtain the path to the ichorCNA conda environment
md5hash = hashlib.md5()
if workflow.conda_prefix:
    conda_prefix = workflow.conda_prefix
else:
    conda_prefix = os.path.abspath(".snakemake/conda")

md5hash.update(conda_prefix.encode())
f = open(CFG["conda_envs"]["ichorcna_run"], 'rb')
md5hash.update(f.read())
f.close()
h = md5hash.hexdigest()

ichorcna_env = conda_prefix + "/" + h + "_" 
ichorcna_dir = ichorcna_env + "/ichorCNA"
ichorcna_scripts_dir = conda_prefix + "/" + h + "_" + "/ichorCNA/scripts/"


rule _ichorcna_install_ichorcna:
    params:
        branch = CFG['options']['ichorcna_branch'],
        directory = ichorcna_dir
    output:
        ichorcna = directory(ichorcna_env + "/lib/R/library/ichorCNA")
    log:
        stdout = CFG["logs"]["inputs"] + "install_ichorCNA.stdout.log",
        stderr = CFG["logs"]["inputs"] + "install_ichorCNA.stderr.log",
    conda:
        CFG['conda_envs']['ichorcna_run']
    shell:
        op.as_one_line("""
        git clone --branch {params.branch} --single-branch https://github.com/broadinstitute/ichorCNA.git {params.directory} &&
        R -q -e "options(timeout=9999999); install.packages('{params.directory}', repos = NULL, type = 'source')"
        """)


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _ichorcna_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam.bai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)


rule _ichorcna_input_pon:
    input:
        pon = CFG["inputs"]["ichorcna_pon"]
    output:
<<<<<<< HEAD
        pon = CFG["dirs"]["inputs"] + "pon/{seq_type}--{genome_build}/bin{binSize}kb/ichorcna_pon.rds"
=======
        pon = CFG["dirs"]["inputs"] + "pon/{seq_type}--{genome_build}/bin{binSize}kb/ichorCNA_PON_{genome_build}_bin{binSize}_full_median.rds"
>>>>>>> main
    run:
        op.absolute_symlink(input.pon, output.pon)


# deeptools to get .bw from .bam
rule _ichorcna_bamCoverage:
    input:
        bam = str(rules._ichorcna_input_bam.output.bam)
    output:
<<<<<<< HEAD
        bw = temp(CFG["dirs"]["bamCoverage"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb.bw")
    params:
        binSize = "{binSize}",
=======
        bw = CFG["dirs"]["bamCoverage"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb.bw"
    params:
>>>>>>> main
        opts = CFG["options"]["bamCoverage"]
    conda: 
        CFG["conda_envs"]["bamCoverage"]
    threads: 
        CFG["threads"]["bamCoverage"]
<<<<<<< HEAD
=======
    group:
        "_ichorcna_bigWig_{tumour_id}"
>>>>>>> main
    resources:
        **CFG["resources"]["bamCoverage"]
    wildcard_constraints: 
        binSize = "|".join(["10", "50", "500", "1000"])
    log:
        CFG["logs"]["bamCoverage"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb.log"
    shell:
        """
<<<<<<< HEAD
            bamCoverage -b {input.bam} --binSize {params.binSize}000 {params.opts} -o {output.bw} -p {threads}
=======
            bamCoverage -b {input.bam} --binSize {wildcards.binSize}000 {params.opts} -o {output.bw} -p {threads}
>>>>>>> main
        """


# Converts bigWig to Wig
rule _ichorcna_bigwigToWig:
    input:
        bw = str(rules._ichorcna_bamCoverage.output.bw)
    output:
        wig = temp(CFG["dirs"]["bigwigToWig"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb_{chrom}.wig")
    conda: 
        CFG["conda_envs"]["bigwigToWig"]
    threads: 
        CFG["threads"]["bigwigToWig"]
    resources:
        **CFG["resources"]["bigwigToWig"]
<<<<<<< HEAD
=======
    group:
        "_ichorcna_Wig_{tumour_id}_{binSize}"
>>>>>>> main
    wildcard_constraints: 
        binSize = "|".join(["10", "50", "500", "1000"]),
        chrom = ".+(?<!--fixed)"
    log:
        CFG["logs"]["bigwigToWig"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb_{chrom}.log"
    shell:
        """
            bigWigToWig {input.bw} {output.wig} -chrom={wildcards.chrom} 
        """


# This function will reformat the wig file to one that can be used for ichorCNA
rule _ichorcna_convert_wig:
    input:
        wig = str(rules._ichorcna_bigwigToWig.output.wig)
    output:
        wig = temp(CFG["dirs"]["convert_wig"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb.{chrom}--fixed.wig")
    conda: 
        CFG["conda_envs"]["convert_wig"]
    threads: 
        CFG["threads"]["convert_wig"]
    resources:
        **CFG["resources"]["convert_wig"]
<<<<<<< HEAD
=======
    group:
        "_ichorcna_Wig_{tumour_id}_{binSize}"
>>>>>>> main
    wildcard_constraints:
        binSize = "|".join(["10", "50", "500", "1000"]),
        chrom = ".+(?<!--fixed)"
    log:
        CFG["logs"]["convert_wig"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb.{chrom}--fixed.log"
    shell:
        """
            echo -e "fixedStep chrom={wildcards.chrom} start=1 step={wildcards.binSize}000 span={wildcards.binSize}000 " > {output.wig} &&
            intersectBed -a <( bedops --chop {wildcards.binSize}000 --header  {input.wig} ) -b {input.wig} -wa -wb | awk '{{print $7}}' >> {output.wig} 2>> {log}
        """


# This function is used to get the wigs of the main chromosomes, which will be stitched together
def get_chrom_wigs(wildcards):
    CFG = config["lcr-modules"]["ichorcna"]
    chrs = reference_files("genomes/" + wildcards.genome_build + "/genome_fasta/main_chromosomes_withY.txt")
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    wigs = expand(
        CFG["dirs"]["convert_wig"] + "{{seq_type}}--{{genome_build}}/bin{{binSize}}kb/{{tumour_id}}.bin{{binSize}}kb.{chrom}--fixed.wig", 
        chrom = chrs)
    return(wigs)


rule _ichorcna_compile_wigs:
    input:
        wigs = get_chrom_wigs
    output:
<<<<<<< HEAD
        wig = temp(CFG["dirs"]["compile_wigs"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb.wig")
=======
        wig = CFG["dirs"]["compile_wigs"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.bin{binSize}kb.wig"
    resources:
        **CFG["resources"]["compile_wigs"]
    group:
        "_ichorcna_run_{tumour_id}"
>>>>>>> main
    shell:
        """
            cat {input.wigs} > {output.wig}
        """


def get_gcwig(wildcards):
    genome_ucsc = op.switch_on_wildcard("genome_build", config["lcr-modules"]["ichorcna"]["options"]["ichorcna_run"]["genome_ucsc"])(wildcards)
    return ichorcna_dir + "/inst/extdata/gc_" + genome_ucsc + "_" + wildcards.binSize + "kb.wig"

def get_mapwig(wildcards):
    genome_ucsc = op.switch_on_wildcard("genome_build", config["lcr-modules"]["ichorcna"]["options"]["ichorcna_run"]["genome_ucsc"])(wildcards)
    return ichorcna_dir + "/inst/extdata/map_" + genome_ucsc + "_" + wildcards.binSize + "kb.wig"

def get_centromere(wildcards):
    centromere = op.switch_on_wildcard("genome_build", config["lcr-modules"]["ichorcna"]["options"]["ichorcna_run"]["centromere"])(wildcards)
    return ichorcna_dir + "/inst/extdata/" + centromere

rule _ichorcna_ichorcna_run:
    input:
        wig = str(rules._ichorcna_compile_wigs.output.wig),
        pon = str(rules._ichorcna_input_pon.output.pon),
        ichorCNA = str(rules._ichorcna_install_ichorcna.output.ichorcna)
    output:
        corrDepth = CFG["dirs"]["ichorcna_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}.correctedDepth.txt",
        param = CFG["dirs"]["ichorcna_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}.params.txt",
        cna = CFG["dirs"]["ichorcna_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}.cna.seg",
        segTxt = CFG["dirs"]["ichorcna_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}.seg.txt",
        seg = CFG["dirs"]["ichorcna_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}.seg",
        plot = CFG["dirs"]["ichorcna_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}/{tumour_id}/{tumour_id}_genomeWide.pdf",
    params:
        ichorcna_scripts_dir = ichorcna_scripts_dir,
        outDir = CFG["dirs"]["ichorcna_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}",
        name = "{tumour_id}",
        gcwig = get_gcwig,
        mapwig = get_mapwig,
        centromere = get_centromere,
        exome = CFG["options"]["ichorcna_run"]["exome"],
        genomeStyle = op.switch_on_wildcard("genome_build", CFG["options"]["ichorcna_run"]["genomeStyle"]),
<<<<<<< HEAD
=======
        genomeBuild = "{genome_build}",
>>>>>>> main
        opts = CFG["options"]["ichorcna_run"]["opts"]
    conda: 
        CFG["conda_envs"]["ichorcna_run"]
    threads: 
        CFG["threads"]["ichorcna_run"]
    resources:
        **CFG["resources"]["ichorcna_run"]
<<<<<<< HEAD
=======
    group:
        "_ichorcna_run_{tumour_id}"
>>>>>>> main
    wildcard_constraints: 
        binSize = "|".join(["10", "50", "500", "1000"])
    log:
        stdout = CFG["logs"]["ichorcna_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.stdout.log",
        stderr = CFG["logs"]["ichorcna_run"] + "{seq_type}--{genome_build}/bin{binSize}kb/{tumour_id}.stderr.log"
    shell:
         op.as_one_line("""
            Rscript {params.ichorcna_scripts_dir}/runIchorCNA.R 
            --id {params.name} 
            --WIG {input.wig} 
            --gcWig {params.gcwig}
            --mapWig {params.mapwig} 
            --normalPanel {input.pon}
            --exons.bed {params.exome}
            --genomeStyle {params.genomeStyle} 
<<<<<<< HEAD
=======
            --genomeBuild {params.genomeBuild}
>>>>>>> main
            --centromere {params.centromere} 
            --outDir {params.outDir} 
            {params.opts}
            > {log.stdout} 2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _ichorcna_output:
    input:
        corrDepth = str(rules._ichorcna_ichorcna_run.output.corrDepth),
        param = str(rules._ichorcna_ichorcna_run.output.param),
        cna = str(rules._ichorcna_ichorcna_run.output.cna),
        seg_txt = str(rules._ichorcna_ichorcna_run.output.segTxt),
        seg = str(rules._ichorcna_ichorcna_run.output.seg),
        plot = str(rules._ichorcna_ichorcna_run.output.plot)
    output:
        corrDepth = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/corrDepth/bin{binSize}kb/{tumour_id}.corrDepth.txt",
        param = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/param/bin{binSize}kb/{tumour_id}.param.txt",
        cna = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/binCNA/bin{binSize}kb/{tumour_id}.cna.seg",
        seg_txt = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/seg_txt/bin{binSize}kb/{tumour_id}.seg.txt",
        seg = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/seg/bin{binSize}kb/{tumour_id}.seg",
        plot = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/plot/bin{binSize}kb/{tumour_id}_genomeWide.pdf"
    run:
        op.relative_symlink(input.corrDepth, output.corrDepth, in_module=True)
        op.relative_symlink(input.param, output.param, in_module=True)
        op.relative_symlink(input.cna, output.cna, in_module=True)
        op.relative_symlink(input.seg_txt, output.seg_txt, in_module=True)
        op.relative_symlink(input.seg, output.seg, in_module=True)
        op.relative_symlink(input.plot, output.plot, in_module=True)


# Generates the target sentinels for each run, which generate the symlinks
rule _ichorcna_all:
    input:
        expand(
            [
                str(rules._ichorcna_output.output.corrDepth),
                str(rules._ichorcna_output.output.param),
                str(rules._ichorcna_output.output.cna),
                str(rules._ichorcna_output.output.seg_txt),
                str(rules._ichorcna_output.output.seg),
<<<<<<< HEAD
                str(rules._ichorcna_output.output.plot)
=======
                str(rules._ichorcna_output.output.plot),
                str(rules._ichorcna_compile_wigs.output.wig)
>>>>>>> main
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            binSize=[CFG["options"]["binsize_kb"]] * len(CFG["runs"]["tumour_sample_id"]))


##### CLEANUP #####

# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
