### Filter VCF files directly based on VAF ###
# requires 3 arguments: input file, output file, and cut-off value for filtering 
# this script will only filter based on VAF, and assumes the VCF is already filtered on PASS value
# currently works with strelka2 outputs only. Support for other variant callers will be added as needed
# Calculation of VAF for both snvs and indels is based on formulas from strelka2 documentation, and strelka2 github issue #3

# example usage: python filter_vcf.py --input <input_file> --output <output_file> --cutoff <value>


#!/usr/bin/env python3

import vcf
import argparse

def main():
    # initiate the parser and handle arguments from command line
    args = parse_args()
    # Open file, this will also read in the header
    vcf_reader = vcf.Reader(open(args.input, 'r'))
    # Open file to write the filterd VCF
    vcf_writer = vcf.Writer(open(args.output, 'w'), vcf_reader)
    cutoff = args.cutoff
    # Determine whenter supplied file is for SNVs or indels and filter accordingly
    print ("Determinig the type of input VCF...")
    if "snv" in str(vcf_reader.metadata['content']):
        print ("VCF file supplied is for snv calls.")
        strelka2_snv(vcf_reader, vcf_writer, cutoff)
    elif "indel" in str(vcf_reader.metadata['content']):
        print ("VCF file supplied is for indel calls.")
        strelka2_indel(vcf_reader, vcf_writer, cutoff)
    else:
        print ("VCF file supplied for filtering is not supported.")
    print ("Done.")
    vcf_writer.close()

# This function filters VCF features based on the VAF for SNV inputs
def strelka2_snv(vcf_reader, vcf_writer, cutoff):
    print ("Extracting VAF and filtering according to the provided cutoff value...")
    for record in vcf_reader:
        if somatic_allele_frequency_snv(record) >= cutoff:
            vcf_writer.write_record(record)
    return vcf_writer

# This function filters VCF features based on the VAF for indel inputs
def strelka2_indel(vcf_reader, vcf_writer, cutoff):
    print ("Extracting VAF and filtering according to the provided cutoff value...")
    for record in vcf_reader:
        if somatic_allele_frequency_indel(record) >= cutoff:
            vcf_writer.write_record(record)
    return vcf_writer

# This function will extract VAF for SNVs from strelka2 outputs
def somatic_allele_frequency_snv(record):
    # Value of FORMAT column $REF + "U"
    reference_allele = str(record.REF).strip("[]")
    refCounts = str(reference_allele+"U")
    # Value of FORMAT column $ALT + "U"
    alternative_allele = str(record.ALT).strip("[]")
    altCounts = str(alternative_allele+"U")
    # First comma-delimited value from $refCounts
    tier1RefCounts = float(record.genotype('TUMOR')[refCounts][0])
    #First comma-delimited value from $altCounts
    tier1AltCounts = float(record.genotype('TUMOR')[altCounts][0])
    # Somatic allele freqeuncy is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)
    allele_freq = tier1AltCounts/(tier1AltCounts+tier1RefCounts)
    return allele_freq

# This function will extract VAF for indels from strelka2 outputs
def somatic_allele_frequency_indel(record):
    # First comma-delimited value from FORMAT/TAR
    tier1RefCounts = float(record.genotype('TUMOR')['TAR'][0])
    # First comma-delimited value from FORMAT/TIR
    tier1AltCounts = float(record.genotype('TUMOR')['TIR'][0])
    # Somatic allele freqeuncy is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)
    allele_freq = tier1AltCounts/(tier1AltCounts+tier1RefCounts)
    return allele_freq


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input",
                        help="Input vcf file to be filtered",
                        required=True)
    parser.add_argument("--output",
                        help="Resulting vcf file after filtering",
                        required=True)
    parser.add_argument("--cutoff",
                        type=float,
                        help="Cutoff value for somatic allele frequency",
                        required=True)

    args, unknown = parser.parse_known_args()

    return args



if __name__ == '__main__':
    main()

