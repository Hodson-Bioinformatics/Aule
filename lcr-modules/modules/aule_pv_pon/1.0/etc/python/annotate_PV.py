'''
Script to annotate found phased variants 
'''

import argparse 
import pysam
import csv
import glob
import re
from Bio.Seq import Seq
import pybedtools
from itertools import combinations
from pyfaidx import Fasta
import pandas as pd

import utilis as ut

# Function to calculate an allele and total depth
def calculate_depth(bam_file, chrom, start, end, ref, alt, min_base_qual, min_mapping_qual):
    
  # Load BAM file
  bam = pysam.AlignmentFile(bam_file, "rb")
  ad = 0
  dp = 0 
  
  start = start-1
  for pileupcolumn in bam.pileup(chrom, int(start), int(end), min_base_quality = int(min_base_qual), min_mapping_quality = int(min_mapping_qual), truncate = True, max_depth = 1000000, ignore_overlaps = True):
    dp = float(pileupcolumn.get_num_aligned())
    seq = pileupcolumn.get_query_sequences()
    ad = seq.count(alt) + seq.count(alt.lower())
  
  bam.close()
  return ad, dp

# Function to annotate all individual PVs with AD and DP
def annotate_depth(bam_file, pv_tab, prefix = "", pv_type = "doublets", min_base_qual = 20, min_mapping_qual = 20):
  
  # Define new columns
  new_columns = ["total_DP", 'AD_pv1', 'DP_pv1', 'AD_pv2', 'DP_pv2']
  new_columns = [prefix + column for column in new_columns]
  
   # Initialize new columns in pv_tab
  for col in new_columns:
    pv_tab[col] = 0

  if pv_type == "triplets":
    extra_columns = ['AD_pv3', 'DP_pv3']
    extra_columns = [prefix + column for column in extra_columns]
    for col in extra_columns:
      pv_tab[col] = 0
  
  # Extract unique genomic positions and alleles
  unique_positions = set()
  for col in ['start_pv1', 'start_pv2'] + (['start_pv3'] if pv_type == "triplets" else []):
    unique_positions.update(zip(pv_tab['seqnames'], pv_tab[col], pv_tab['ref_' + col[-3:]], pv_tab['alt_' + col[-3:]]))

  # Calculate metrics for unique positions
  depth_metrics = {}
  for seqname, position, ref, alt in unique_positions:
    ad, dp = calculate_depth(bam_file, seqname, position, position, ref, alt, min_base_qual, min_mapping_qual)
    depth_metrics[(seqname, position, ref, alt)] = (ad, dp)
      
  # Assign calculated metrics to PVs
  for index, row in pv_tab.iterrows():
    # PV1
    key = (row['seqnames'], row['start_pv1'], row['ref_pv1'], row['alt_pv1'])
    pv_tab.at[index, new_columns[1]], pv_tab.at[index, new_columns[2]] = depth_metrics[key]

    # PV2
    key = (row['seqnames'], row['start_pv2'], row['ref_pv2'], row['alt_pv2'])
    pv_tab.at[index, new_columns[3]], pv_tab.at[index, new_columns[4]] = depth_metrics[key]

    if pv_type == "triplets":
      # PV3
      key = (row['seqnames'], row['start_pv3'], row['ref_pv3'], row['alt_pv3'])
      pv_tab.at[index, extra_columns[0]], pv_tab.at[index, extra_columns[1]] = depth_metrics[key]
      
  return pv_tab


# Check if genomic positions overlap with regions in a BED file.
def check_overlap(chromosomes, positions, bed_file_path):
    
    # Create a BEDTools object from the provided chromosomes and positions
    query_bed = pybedtools.BedTool([(chrom, pos, pos+1) for chrom, pos in zip(chromosomes, positions)])

    # Load the BED file
    bed = pybedtools.BedTool(bed_file_path)

    # Check for overlaps
    overlap_results = query_bed.intersect(bed, u=True)

    # Create a set of tuples for the overlaps
    overlaps = set([(entry[0], int(entry[1])) for entry in overlap_results])

    # Check each position for overlap
    overlap_vector = [1 if (chrom, pos) in overlaps else 0 for chrom, pos in zip(chromosomes, positions)]

    return overlap_vector

# Master function performing annotations for PVS 
def annotate_pv(bam, pv_in, pv_type, panel, pv_out, min_base_qual, min_mapping_qual):
    
  # Load PVs
  pvs = pd.read_csv(pv_in, delimiter='\t')

  # Check if a PVs is on target
  pvs["on_target_pv1"] = check_overlap(pvs['seqnames'], pvs['start_pv1'], panel)
  pvs["on_target_pv2"] = check_overlap(pvs['seqnames'], pvs['start_pv2'], panel)
  if pv_type == "triplets":
    pvs["on_target_pv3"] = check_overlap(pvs['seqnames'], pvs['start_pv3'], panel)
    
  # Filter - at least one PV compnent on target 
  if pv_type == "doublets":
    tumour_pv = pvs[(pvs['on_target_pv1'] == 1) | (pvs['on_target_pv2'] == 1)]
  elif pv_type == "triplets":
    tumour_pv = pvs[(pvs['on_target_pv1'] == 1) | (pvs['on_target_pv2'] == 1) | (pvs['on_target_pv3'] == 1)]
    
  pv_n = len(pvs)
  print(["Number of PV to process:", pv_n])
  
  # Annotate PVs with total depth and AD for individual PV components
  print(["Annotating depth"])
  annotate_depth(bam, pvs, pv_type, min_base_qual, min_mapping_qual)
  
  # Save doublets
  pvs.to_csv(pv_out, sep='\t', index=False)
  print(tumour_pv)
  
  return None


def main():

    # Parse arguments 
    parser = argparse.ArgumentParser()
    parser.add_argument('-bam', metavar='<bam>', dest="bam", help="input bam file")
    parser.add_argument('-input_doublets', metavar='<input_doublets>', dest="input_doublets", help="input doublets")
    parser.add_argument('-input_triplets', metavar='<input_triplets>', dest="input_triplets", help="input triplets")
    parser.add_argument('-output_doublets', metavar='<output_doublets>', dest="output_doublets", help="output doublets")
    parser.add_argument('-output_triplets', metavar='<output_triplets>', dest="output_triplets", help="output triplets")
    parser.add_argument('-panel', metavar='<panel>', dest="panel", help="output panel")
    parser.add_argument('-min_mapping_qual', metavar='<min_mapping_qual>', dest="min_mapping_qual", help="min MAPQ")
    parser.add_argument('-min_base_qual', metavar='<min_base_qual>', dest="min_base_qual", help="min BQ")
    
    args = parser.parse_args()
             
    # Iterate through the file saving variants metrics to the output file
    annotate_pv(args.bam, args.input_doublets, "doublets",  args.panel, args.output_doublets, args.min_base_qual, args.min_mapping_qual)
    annotate_pv(args.bam, args.input_triplets, "triplets",  args.panel, args.output_triplets, args.min_base_qual, args.min_mapping_qual)

if __name__ == "__main__":
	main()
