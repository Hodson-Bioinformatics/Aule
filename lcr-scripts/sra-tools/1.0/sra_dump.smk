import pandas as pd
import os

configfile: "sra_dump.yaml"

# Load the accession numbers from file. 
accessions = pd.read_csv(config["accessions"], names = ["accession"])

# Confirm the output file type is correctly specified. 
out_file_type = config["out_file_type"]
possible_types = ["bam", "cram", "fastq"]
possible_type_string = ", ".join(possible_types)
assert out_file_type in possible_types, (
    f"Unrecognized out_file_type {out_file_type}. \n"
    f"Please specify one of {possible_type_string}. "
)

wildcard_constraints:
    accession = "|".join(accessions["accession"].tolist())

# Create an argument to use the .ngc file if it is specified. 
def get_ngc_arg(ngc = config["params"]["ngc"]):
    if ngc not in ["__UPDATE__", ""]:
        assert os.path.exists(ngc), (
            f"ngc file {ngc} does not exist. \n"
            "Please specify a .ngc file or leave config[\"params\"][\"ngc\"] blank if the repository is public. "
        )
        ngc_param = f"--ngc {ngc}"
    else:
        ngc_param = ""
    return ngc_param

# Prefetch the files. This is done because it is much faster than using sam-dump or fast(er)q-dump directly, and can be resumed if it gets interrupeted. 
rule prefetch:
    output:
        fetched = temp(directory("{accession}"))
    log:
        log = "{accession}.log"
    params:
        ngc = get_ngc_arg(),
        opts = config["params"]["prefetch"]
    conda:
        config["conda"]["sra-tools"]
    threads: config["threads"]["prefetch"]
    resources:
        **config["resources"]["prefetch"]
    group: "sra-dump"
    shell:
        "prefetch -X {resources.disk_mb} {params.opts} {params.ngc} {wildcards.accession} > {log.log} 2>&1"

# Convert the prefetched file to fastq if desired output is fastq. 
rule fastq_dump:
    input:
        fetched = str(rules.prefetch.output.fetched)
    output:
        fq1 = temp("{accession}_1.fastq"),
        fq2 = temp("{accession}_2.fastq"),
        complete = touch("{accession}.fastq.complete")
    params:
        opts = config["params"]["fastq_dump"]
    conda:
        config["conda"]["sra-tools"]
    threads: config["threads"]["fastq_dump"]
    resources:
        **config["resources"]["fastq_dump"]
    group: "sra-dump"
    shell:
        "fasterq-dump -e {threads} {params.opts} {wildcards.accession}"


# Compress the output fastq file. 
rule gzip_fastq:
    input:
        fastq = "{accession}_{read}.fastq",
        complete = "{accession}.fastq.complete"
    output:
        fastq = "{accession}_{read}.fastq.gz", 
        complete = touch("{accession}_{read}.gzipped")
    threads: 1
    resources:
        **config["resources"]["gzip_fastq"]
    wildcard_constraints:
        read = "1|2"
    group: "sra-dump"
    shell:
        "gzip -c {input.fastq} > {output.fastq}"



# Obtain the reference genome argument for cram compression. 
def get_genome_arg(
    out_type = config["out_file_type"],
    ref_fasta = config["params"]["ref"]
):
    if out_type == "cram":
        assert os.path.exists(ref_fasta), (f"Reference fasta file {ref_fasta} does not exist. ")
        param = f"-CS -T {ref_fasta}"
    elif out_type == "bam":
        param = "-bS"
    else:
        param = ""
    return(param)


# Convert prefetched file to sam if specified. 
rule sam_dump:
    input:
        fetched = str(rules.prefetch.output.fetched)
    output:
        sam = pipe("{accession}.sam")
    params:
        ngc = get_ngc_arg(),
        opts = config["params"]["sam-dump"]
    conda:
        config["conda"]["sra-tools"]
    threads: 1
    resources:
        mem_mb = 10000
    group: "sra-dump"
    shell:
        "sam-dump {params.opts} {wildcards.accession} > {output.sam}"

# Compress sam to bam or cram. 
rule sam_to_bam:
    input:
        sam = str(rules.sam_dump.output.sam)
    output:
        bam = "{accession}.hg19." + config["out_file_type"],
        complete = touch("{accession}.complete")
    params:
        bam_cram = get_genome_arg(),
        opts = config["params"]["samtools"]
    conda:
        config["conda"]["sra-tools"]
    threads: config["threads"]["sam_to_bam"]
    resources:
        **config["resources"]["sam_to_bam"]
    group: "sra-dump"
    shell:
        "samtools view {params.bam_cram} -@ {threads} {input.sam} > {output.bam}"


# Generate an index file. 
rule index_bam:
    input:
        bam = str(rules.sam_to_bam.output.bam)
    output:
        complete = touch("{accession}.indexed")
    conda:
        config["conda"]["sra-tools"]
    threads: config["threads"]["index_bam"]
    resources:
        **config["resources"]["index_bam"]
    group: "sra-dump"
    shell:
        "samtools index -@ {threads} {input.bam}"


rule all:
    input:
        expand(
            [
                str(rules.sam_to_bam.output.bam),
                str(rules.sam_to_bam.output.complete),
                str(rules.index_bam.output.complete)
            ],
            accession = accessions["accession"]
        ) if out_file_type in ["cram", "bam"] else expand(
            expand[
                str(rules.gzip_fastq.output.fastq), 
                str(rules.gzip_fastq.output.complete)
            ],
            accession = accessions["accession"],
            read = ["1", "2"]
        )