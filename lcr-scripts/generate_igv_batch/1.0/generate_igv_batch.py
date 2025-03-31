#!/usr/bin/env python3

import re
import os
import glob
import argparse

import vcf
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(prog='generate_igv_batch')
parser.add_argument('--bam', nargs='*',  help='list of paths to bam files you want to include in the snapshot (e.g. tumour, normal, RNA-seq)',required=True)
parser.add_argument('--input', nargs='?', help='path to vcf file, either uncompressed or compressed, bed/bedpe, or maf')
parser.add_argument('--mode', nargs='?', help='Specify whether snapshots will be of one region ("single") or two regions ("pair" e.g. for SVs). This only works for vcf or bed/bedpe.',required=True)
parser.add_argument('--batchfile',nargs='?',help="name of igv batch file to create",required=True)
parser.add_argument('--basename',nargs='?',help="name of sample or patient for base filename for each snapshot",required=True)
parser.add_argument('--snapshot_dir',nargs='?',help="where all the snapshots will be created (this directory must already exist)",required=True)
parser.add_argument('--pad',nargs='?',default=200,type=int,help="how many nucleotides padding to add around position")
parser.add_argument('--genome_build',help="specify whether to use hg19 or hg38",required=True)
parser.add_argument('--max_height',default=400,type=int,help="maximum panel height in IGV")
parser.add_argument('--pass_only',action='store_true',help="flag to force snapshots for only the lines with FILTER=PASS")
parser.add_argument('--geneList', nargs = '?', help="add a list of genes to examine")

args = parser.parse_args()

print(args)

#detects if input was submitted
if not args.input:
    print("must supply either a bed/bedpe, maf, or vcf file as input")
    exit()

#converts chr, start, end in different file types into a format that is readable by IGV (i.e. chr:start-end)
def get_regions(region_mode=args.mode):
    if args.input.endswith(('.vcf', '.vcf.gz')): #for vcf input
        vcf_reader = vcf.Reader(filename=args.input,compressed=True)
        regions = []
        for record in vcf_reader:
            if args.pass_only:
                if len(record.FILTER)>0:
                    continue
            region1 = "{}:{}-{}".format(record.CHROM, record.POS-args.pad,record.POS+args.pad)
            if region_mode == "single":
                regions.append((region1,))
            elif region_mode == "pair":
                region2 = "{}:{}-{}"
                if 'END' in record.INFO:
                    region2=region2.format(record.CHROM,record.INFO['END']-args.pad,record.INFO['END']+args.pad)
                else:
                    region2=region2.format(record.ALT[0].chr,record.ALT[0].pos-args.pad,record.ALT[0].pos+args.pad)
                regions.append((region1,region2))
        return regions
    elif args.input.endswith(('.bed', '.bedpe')): #for bed input
        regions = []
        with open(args.input) as bed_file_raw:
            if region_mode == "single":
                for line in bed_file_raw: #takes line by line
                    fields = line.rstrip("\n").split("\t")
                    chrom, start, end = fields[:3] #only col 1-3 needed
                    output = f"{chrom}:{start}-{end}" #output is the formatted genomic coordinate
                    regions.append((output,))
            elif region_mode == "pair":
                for line in bed_file_raw: #takes line by line
                    fields = line.rstrip("\n").split("\t")
                    chrom1, start1, end1, chrom2, start2, end2 = fields[:6] #only col 1-6 needed
                    output1 = f"{chrom1}:{start1}-{end1}"
                    output2 = f"{chrom2}:{start2}-{end2}"
                    regions.append((output1,output2))
        return regions
    elif args.input.endswith(('.maf')): #for maf input
        regions = [] # creates empty variable
        maf_file_raw = pd.read_csv(args.input, skiprows=1, sep='\t') # loads file as a data.frame #skips line 1
        maf_file_raw['coord'] = maf_file_raw['Chromosome'].astype(str)+':'+maf_file_raw['Start_Position'].astype(str)+'-'+maf_file_raw['End_Position'].astype(str) # creates a new column of formatted coordinates
        if not args.geneList:
            print("No gene list. Continuing...")
        else:
            genes = pd.read_csv(args.geneList, header = None) # opens txt; no header
            genes = genes[0].tolist() # turns dataframe to genelist
            maf_file_raw = maf_file_raw[maf_file_raw['Hugo_Symbol'].isin(genes)] #filters maf for genes of interest
        maf_file = maf_file_raw[['coord','Hugo_Symbol']]
        regions = np.array(maf_file)
        return regions
regions = get_regions(args.mode)

### generates the formatted batch file
# batch commands #  not written online
# setLogScale(true | false)
# remove TrackName
outfile = open(args.batchfile,"w")
for region_list in regions:
    for bamfile in args.bam:
        outfile.write("load {}\n".format(bamfile))
    outfile.write("maxPanelHeight {}\n".format(args.max_height))
    outfile.write("snapshotDirectory {}\n".format(args.snapshot_dir))
    outfile.write("genome {}\n".format(args.genome_build))
    outfile.write("goto ")
    if args.input.endswith(('.maf')):
        for region in region_list[[0]]:
            outfile.write("{}".format(region))
    else:
        for region in region_list:
            outfile.write("{}".format(region))
    outfile.write("\n")
    outfile.write("sort\n")
    outfile.write("collapse\n")
    for bamfile in args.bam:
        sample = os.path.basename(bamfile)
        outfile.write("remove '{} Junctions'\n".format(sample))
    if len(region_list) >1:
        outfile.write("snapshot {}-{}-{}.png\n".format(args.basename,region_list[0],region_list[1]))
    else:
        outfile.write("snapshot {}-{}.png\n".format(args.basename,region_list[0]))
    outfile.write("new\n")
outfile.write("exit\n")




