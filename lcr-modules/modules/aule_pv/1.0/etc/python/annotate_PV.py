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
def calculate_depth(bam_file, chrom, start, end, ref, alt):
    
  # Load BAM file
  bam = pysam.AlignmentFile(bam_file, "rb")
  ad = 0
  dp = 0 
  
  start = start-1
  for pileupcolumn in bam.pileup(chrom, start, end, min_base_quality = 20, min_mapping_quality = 20, truncate = True, max_depth = 1000000, ignore_overlaps = True):
    dp = float(pileupcolumn.get_num_aligned())
    seq = pileupcolumn.get_query_sequences()
    ad = seq.count(alt) + seq.count(alt.lower())
  
  bam.close()
  return ad, dp

def calculate_total_depth_pairs(bam_file, chrom, start, end):
    
  # Load BAM file
  bam = pysam.AlignmentFile(bam_file, "rb")
  dp = 0 
  
  start = start-1
  for pileupcolumn in bam.pileup(chrom, start, end, min_base_quality = 20, min_mapping_quality = 20, truncate = True, max_depth = 1000000, ignore_overlaps = True):
    dp = float(pileupcolumn.get_num_aligned())
                    
  bam.close()
  return dp

# Function to annotate all individual PVs with AD and DP
def annotate_depth(bam_file, pv_tab, prefix = "", pv_type = "doublets"):
  
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
    ad, dp = calculate_depth(bam_file, seqname, position, position, ref, alt)
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

    ## Total depth (measured by the number of read pairs)
    #pv_tab.at[index, new_columns[0]] = calculate_total_depth_pairs(bam_file, row['seqnames'], row['end5'], row['end3'])

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

# Extract 5nt context around PVs
def extract_5nt_context(chromosomes, positions, genome):
  context_sequences = []
  
  for chrom, pos in zip(chromosomes, positions):
    start = pos - 3  # Adjust to get two bases before the position
    end = pos + 2    # Adjust to get two bases after the position
    sequence = genome[chrom][start:end].seq
    context_sequences.append(sequence)

  return context_sequences
  
  
# Annotate with PON
def extract_fraction_from_vcf(chromosomes, positions, vcf_file_path, fraction_info):
  
  # Load VCF file
  vcf = pysam.VariantFile(vcf_file_path)
  fraction_values = []
  
  for chrom, pos in zip(chromosomes, positions):
      fraction = 0  # Default value if no overlap
      for rec in vcf.fetch(chrom, pos-1, pos):
          if fraction_info in rec.info:
            fraction = rec.info[fraction_info]
            if not isinstance(fraction, float):
              fraction = fraction[0]
            break # Stop after the first matching record
      fraction_values.append(fraction)
  vcf.close()
  return fraction_values


# Master function performing annotations for PVS 
def annotate_pv(tumour_bam, normal_bam, tumour_pv_file, normal_pv_file, pv_type, fasta, pv_pon_file, pon_snv, gnomad, panel, pv_out):
    
  # Load PVs
  tumour_pv = pd.read_csv(tumour_pv_file, delimiter='\t')
  normal_pv = pd.read_csv(normal_pv_file, delimiter='\t')
  
  # Load PV PON
  if pv_type == "doublets":
    pv_pon = pd.read_csv(pv_pon_file, delimiter='\t')
  elif pv_type == "triplets":
    pv_pon = pd.read_csv(pv_pon_file, delimiter='\t')
  
  # Check if a PVs is on target
  tumour_pv["on_target_pv1"] = check_overlap(tumour_pv['seqnames'], tumour_pv['start_pv1'], panel)
  tumour_pv["on_target_pv2"] = check_overlap(tumour_pv['seqnames'], tumour_pv['start_pv2'], panel)
  if pv_type == "triplets":
    tumour_pv["on_target_pv3"] = check_overlap(tumour_pv['seqnames'], tumour_pv['start_pv3'], panel)
    
  # Filter - at least one PV compnent on target 
  if pv_type == "doublets":
    tumour_pv = tumour_pv[(tumour_pv['on_target_pv1'] == 1) | (tumour_pv['on_target_pv2'] == 1)]
  elif pv_type == "triplets":
    tumour_pv = tumour_pv[(tumour_pv['on_target_pv1'] == 1) | (tumour_pv['on_target_pv2'] == 1) | (tumour_pv['on_target_pv3'] == 1)]
    
  pv_n = len(tumour_pv)
  print(["Number of PV to process:", pv_n])
  
  ## DOUBLETS
  # Annotate PVs with total depth and AD for individual PV components
  print(["Annotating depth"])
  annotate_depth(tumour_bam, tumour_pv, pv_type = pv_type)
  annotate_depth(normal_bam, tumour_pv, pv_type = pv_type, prefix = "germline_")
  
  # Load genome 
  genome = Fasta(fasta)
  
  # Annotate variants with 5 nuclotides context 
  tumour_pv["context_pv1"] = extract_5nt_context(tumour_pv['seqnames'], tumour_pv['start_pv1'], genome)
  tumour_pv["context_pv2"] = extract_5nt_context(tumour_pv['seqnames'], tumour_pv['start_pv2'], genome)
  
  # Annotate with PON
  tumour_pv["pon_pv1"] = extract_fraction_from_vcf(tumour_pv['seqnames'], tumour_pv['start_pv1'], pon_snv, fraction_info = "FRACTION")
  tumour_pv["pon_pv2"] = extract_fraction_from_vcf(tumour_pv['seqnames'], tumour_pv['start_pv2'], pon_snv, fraction_info = "FRACTION")
  
  # Annotate with GNOMAD
  tumour_pv["gnomad_pv1"] = extract_fraction_from_vcf(tumour_pv['seqnames'], tumour_pv['start_pv1'], gnomad, fraction_info = "AF")
  tumour_pv["gnomad_pv2"] = extract_fraction_from_vcf(tumour_pv['seqnames'], tumour_pv['start_pv2'], gnomad, fraction_info = "AF")
  
  # Annotate with PV PON 
  tumour_pv = pd.merge(tumour_pv, normal_pv[['PV_ID', 'AD_qual', 'AD_total']], on='PV_ID', how='left', suffixes=('_tumour', '_germline'))
  tumour_pv['AD_qual_germline'].fillna(0, inplace=True)
  tumour_pv['AD_total_germline'].fillna(0, inplace=True)
  
  # Operations for triplets
  if pv_type == "triplets":
    tumour_pv["context_pv3"] = extract_5nt_context(tumour_pv['seqnames'], tumour_pv['start_pv3'], genome)
    tumour_pv["pon_pv3"] = extract_fraction_from_vcf(tumour_pv['seqnames'], tumour_pv['start_pv3'], pon_snv, fraction_info = "FRACTION")
    tumour_pv["gnomad_pv3"] = extract_fraction_from_vcf(tumour_pv['seqnames'], tumour_pv['start_pv3'], gnomad, fraction_info = "AF")
    
  # Annotate tumour doublets with normal AD
  pv_pon.rename(columns={'meanPVAF': 'PON_PV_meanPVAF', 'PVAF_total': 'PON_PV_PVAFtotal', 'AD_average': 'PON_PV_meanAD', 'N':'PON_PV_Nsamples'}, inplace=True)
  tumour_pv = pd.merge(tumour_pv, pv_pon[['PV_ID', 'PON_PV_meanPVAF', 'PON_PV_PVAFtotal', 'PON_PV_meanAD', 'PON_PV_Nsamples']], on='PV_ID', how='left')
  tumour_pv['PON_PV_meanPVAF'].fillna(0, inplace=True)
  tumour_pv['PON_PV_PVAFtotal'].fillna(0, inplace=True)
  tumour_pv['PON_PV_meanAD'].fillna(0, inplace=True)
  tumour_pv['PON_PV_Nsamples'].fillna(0, inplace=True)
  
  # Save doublets
  tumour_pv.to_csv(pv_out, sep='\t', index=False)
  print(tumour_pv)
  
  return None


def main():

    # Parse arguments 
    parser = argparse.ArgumentParser()
    parser.add_argument('-tumour_bam', metavar='<tumour_bam>', dest="tumour_bam", help="input tumour bam file")
    parser.add_argument('-normal_bam', metavar='<normal_bam>', dest="normal_bam", help="input normal bam file")
    parser.add_argument('-tumour_pv', metavar='<tumour_pv>', dest="tumour_pv", help="PV in tumour sample")
    parser.add_argument('-normal_pv', metavar='<normal_pv>', dest="normal_pv", help="PV in normal sample")
    parser.add_argument('-pv_type', metavar='<pv_type>', dest="pv_type", help="PV type [doublets|triplets]")
    parser.add_argument('-fasta', metavar='<fasta>', dest="fasta", help="reference fasta")
    parser.add_argument('-pon_pv', metavar='<pon_pv>', dest="pon_pv", help="vcf with panel of normals (phased variants)")
    parser.add_argument('-pon_snv', metavar='<pon_snv>', dest="pon_snv", help="vcf with panel of normals (snvs)")
    parser.add_argument('-gnomad', metavar='<gnomad>', dest="gnomad", help="vcf with population SNPs")
    parser.add_argument('-panel', metavar='<panel>', dest="panel", help="BED file with panel coordinates")
    parser.add_argument('-pv_out', metavar='<pv_out>', dest="pv_out", help="annotated phased variants")

    args = parser.parse_args()
               
    # Iterate through the file saving variants metrics to the output file
    annotate_pv(args.tumour_bam, args.normal_bam, args.tumour_pv, args.normal_pv, args.pv_type, args.fasta, args.pon_pv, args.pon_snv, args.gnomad, args.panel, args.pv_out)

if __name__ == "__main__":
	main()
