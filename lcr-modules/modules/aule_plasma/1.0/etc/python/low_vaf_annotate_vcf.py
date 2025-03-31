'''
Script to annotate a VCF file with features for ultra-sensitive variant calling. 
It is designed to be an extension to a SLMS3 workflow, but any vcf file can be processed
Each variant in the provided VCF file will be annotated wit the following features:
  - Allele depth 
'''

import argparse 
import pysam
import csv
import glob
import re
import numpy as np
import pandas as pd
import logging
from itertools import combinations


# Get variant classification
def get_variant_class(ref, alts):

  var_types = []
  var_length = []
  for alt in alts:
    if len(alt) == 1 and len(ref) > 1:
      var_types.append("DEL")
      var_length.append(1)
    elif len(alt) == 1 and len(ref) == 1:
      var_types.append("SNP")
      var_length.append(1)
    elif len(alt) > 1 and len(ref) == 1:
      var_types.append("INS")
      var_length.append(len(alt))
    elif len(alt) == 2 and len(ref) == 2:
      var_types.append("DNP")
      var_length.append(2)
    elif len(alt) > 2 and len(ref) > 2:
      var_types.append("MNP")
      var_length.append(int(len(alt)))
      
  return({"var_types": var_types, "var_length": var_length})


# Parse basic variant metrics 
def parse_variant(variant):
  
  # Define coordinates 
  chrom = variant.chrom
  ref = variant.ref
  alts = list(variant.alts)
  start = variant.pos-1 # The actual start (pysam uses 0-coordinate system)
  var_metrics = get_variant_class(ref, alts)
  end = [x + start for x in var_metrics["var_length"]]
  
  return({"chrom": chrom, 
          "start": start, 
          "end": end, 
          "ref": ref, 
          "alts": alts, 
          "var_types": var_metrics["var_types"], 
          "var_length": var_metrics["var_length"]})
          
# Compute mean values per base 
def compute_mean_per_base(sum_values, base_counts):
  
  mean_values = sum_values
  for key in sum_values:
    try:
        mean_values[key] = sum_values[key] / base_counts[key]
    except ZeroDivisionError:
        mean_values[key] = None  # or you can choose another way to handle it, like setting it to 0
  
  return(mean_values)

# Collect metrics at the base location in a given read
def collect_position_metrics(pileupread, pos):
  
  #print(pileupread.alignment.query_qualities[pos])
  this_qual = pileupread.alignment.query_qualities[pos] # Get the base quality
  this_family_size = pileupread.alignment.get_tag("cd")[pos] # Get the family size at this location
  this_umi_discordance = pileupread.alignment.get_tag("ce")[pos] # Is the location UMI discordant?
  this_fragmlen = abs(pileupread.alignment.template_length) # Get the fragment length
    
  return({"this_qual": this_qual, "this_family_size": this_family_size, "this_umi_discordance": this_umi_discordance,"this_fragmlen": this_fragmlen })

  
# Repalce None with 0  
def replace_none_with_zero(dictionary):
    return {k: 0 if v is None else v for k, v in dictionary.items()}
  

# Sumamrise pileup 
def summarise_pileup(pileupcolumn, var_type, var_alt, var_ref):
  
  alt_count = 0
  ref_count = 0
  base_counts = {"ALT": 0 ,"REF": 0, "N": 0, "A": 0, "T": 0, "G": 0, "C": 0, "INS": 0, "DEL": 0, "MNP": 0, "VAR": 0}
  nontrivial_families = {"ALT": 0 ,"REF": 0 , "N": 0, "A": 0, "T": 0, "G": 0, "C": 0, "INS": 0, "DEL": 0, "MNP": 0, "VAR": 0}
  reverse_reads = {"ALT": 0 ,"REF": 0 , "N": 0, "A": 0, "T": 0, "G": 0, "C": 0, "INS": 0, "DEL": 0, "MNP": 0, "VAR": 0}
  umi_discordant = {"ALT": 0 ,"REF": 0 , "N": 0, "A": 0, "T": 0, "G": 0, "C": 0, "INS": 0, "DEL": 0, "MNP": 0, "VAR": 0}
  max_familysize = {"ALT": 0 ,"REF": 0 , "N": 0, "A": 0, "T": 0, "G": 0, "C": 0, "INS": 0, "DEL": 0, "MNP": 0, "VAR": 0}
  mean_bq = {"ALT": 0 , "REF": 0 ,"N": 0, "A": 0, "T": 0, "G": 0, "C": 0, "INS": 0, "DEL": 0, "MNP": 0, "VAR": 0}
  mean_tlen = {"ALT": 0 , "REF": 0 , "N": 0, "A": 0, "T": 0, "G": 0, "C": 0, "INS": 0, "DEL": 0, "MNP": 0, "VAR": 0}
  
  all_read_names = []
  
  n = 0
  pileup_seq = pileupcolumn.get_query_sequences(add_indels=True)
  pileup_qual = pileupcolumn.get_query_qualities()

  
  # Iterate through reads
  for pileupread in pileupcolumn.pileups:
    
    this_read_name = pileupread.alignment.query_name
    if not this_read_name in all_read_names:
      all_read_names.append(this_read_name)
      
      pos = pileupread.query_position_or_next # Position in the read
      #print(pos, pileupread.alignment.query_alignment_length, len(pileupread.alignment.query_qualities))
      #print(pileupread.alignment.query_qualities[pos])
      ### SNPs
      if var_type == "SNP":
        if not pileupread.is_refskip: # query position is None if is_del or is_refskip is set.
          
          # Identify base sequence at this position in this read
          if "+" in pileup_seq[n]: # Identify reads with insertions at this position
            base_counts["INS"] +=1 
            this_base = "INS"
          elif "-" in pileup_seq[n]:
            base_counts["DEL"] +=1 
            this_base = "DEL"
          elif len(var_alt) > 1:
            base_counts["MNP"] +=1 
            this_base = "MNP"
          else:
            this_base = pileupread.alignment.query_sequence[pos].upper()
            base_counts[this_base] +=1
            if this_base == var_alt:
              alt_count +=1
    
      ### INSERTIONS
      elif var_type == "INS":
        if not pileupread.is_refskip:  
          
          # Identify base sequence at this position in this read
          if "+" in pileup_seq[n]: # Identify reads with insertions at this position
            this_base = pileupread.alignment.query_sequence[pos:pos+len(var_alt)].upper()
            if this_base == var_alt: # Count insertion at insertion site only if fully matched with alt
              base_counts["INS"] +=1 
              alt_count +=1
              this_base = "INS"
            else:
              this_base = "VAR"
          elif "-" in pileup_seq[n]:
            base_counts["DEL"] +=1 
            this_base = "DEL"
          elif len(var_alt) > 1:
            base_counts["MNP"] +=1 
            this_base = "MNP"
          else:
            this_base = pileupread.alignment.query_sequence[pos].upper()
            base_counts[this_base] +=1 
      
      ### DELETIONS
      elif var_type == "DEL":
        del_length = len(var_ref)-1
        # Identify base sequence at this position in this read
       
        if "-" in pileup_seq[n]:
          del_bases = pileup_seq[n].count('N')
          if del_bases == del_length: # Count insertion at insertion site only if fully matched with alt
            base_counts["DEL"] +=1 
            this_base = "DEL"
            alt_count +=1
          else:
            this_base = "VAR"  
        elif "+" in pileup_seq[n]:
          base_counts["INS"] +=1 
          this_base = "INS"
        elif len(var_alt) > 1:
            base_counts["MNP"] +=1 
            this_base = "MNP"
        else:
          this_base = pileupread.alignment.query_sequence[pos].upper()
          base_counts[this_base] +=1
      
      ### MUTLIPLE BASE SUBSTITUTIONS    
      elif var_type == "DNP" or var_type == "MNP":
        if len(var_alt) == 1:
          this_base = pileupread.alignment.query_sequence[pos].upper()
          base_counts[this_base] +=1
        elif "-" in pileup_seq[n]:
          base_counts["DEL"] +=1 
          this_base = "DEL"
        elif "+" in pileup_seq[n]:
          base_counts["INS"] +=1 
          this_base = "INS"
        else:
          this_base = pileupread.alignment.query_sequence[pos:pos+len(var_alt)].upper()
          if this_base == var_alt: # Count only if fully matched with alt
            alt_count +=1
            this_base = "MNP"
          else:
            this_base = "VAR"
        
      # Collect base metrics
      base_metrics = collect_position_metrics(pileupread, pos)
      
      # Save metrics for ALT and REF, if aplicable 
      if this_base == var_ref[0]:
        ref_count +=1
        mean_bq["REF"] += base_metrics["this_qual"]
        mean_tlen["REF"] += base_metrics["this_fragmlen"]
        max_familysize["REF"] = max(max_familysize["REF"], base_metrics["this_family_size"])
        if base_metrics["this_family_size"] > 1:
          nontrivial_families["REF"] += 1        # Count nontrivial famili sizes
        if base_metrics["this_umi_discordance"] > 0: # Count discordant umi families
          umi_discordant["REF"] += 1
        if pileupread.alignment.is_reverse:          # Count reverse reads
          reverse_reads["REF"] += 1   
      elif this_base == var_alt:
        mean_bq["ALT"] += base_metrics["this_qual"]
        mean_tlen["ALT"] += base_metrics["this_fragmlen"]
        max_familysize["ALT"] = max(max_familysize["ALT"], base_metrics["this_family_size"])
        if base_metrics["this_family_size"] > 1:
          nontrivial_families["ALT"] += 1        # Count nontrivial famili sizes
        if base_metrics["this_umi_discordance"] > 0: # Count discordant umi families
          umi_discordant["ALT"] += 1
        if pileupread.alignment.is_reverse:          # Count reverse reads
          reverse_reads["ALT"] += 1   

      # Save per base metrics
      mean_bq[this_base] += base_metrics["this_qual"]
      mean_tlen[this_base] += base_metrics["this_fragmlen"]
      max_familysize[this_base] = max(max_familysize[this_base], base_metrics["this_family_size"])
      if base_metrics["this_family_size"] > 1:
        nontrivial_families[this_base] += 1        # Count nontrivial famili sizes
      if base_metrics["this_umi_discordance"] > 0: # Count discordant umi families
        umi_discordant[this_base] += 1
      if pileupread.alignment.is_reverse:          # Count reverse reads
        reverse_reads[this_base] += 1      
    
      # Keep track of reads 
      n += 1
  
  # Compute mean values
  base_counts["REF"] = ref_count
  base_counts["ALT"] = alt_count
  
  mean_bq = compute_mean_per_base(mean_bq, base_counts)
  mean_tlen = compute_mean_per_base(mean_tlen, base_counts)
  
  # Compute vaf
  vaf = alt_count/n
  
  # Replace None with 0 
  nontrivial_families = replace_none_with_zero(nontrivial_families)
  reverse_reads = replace_none_with_zero(reverse_reads)
  umi_discordant = replace_none_with_zero(umi_discordant)
  max_familysize = replace_none_with_zero(max_familysize)
  mean_bq = replace_none_with_zero(mean_bq)
  mean_tlen = replace_none_with_zero(mean_tlen)
  
  return({"alt_counts": alt_count,
          "ref_counts": ref_count,
          "total_counts": n,
          "vaf": vaf, 
          "base_counts":base_counts, 
          "nontrivial_families": nontrivial_families, 
          "reverse_reads": reverse_reads, 
          "umi_discordant": umi_discordant, 
          "max_familysize": max_familysize, 
          "mean_bq": mean_bq, 
          "mean_tlen": mean_tlen})
  

# Extract variant metrics from the bam files
def get_variant_metrics(variant, bam, min_base_qual, min_mapping_qual, new_header):
  
  variants_metrics = []
  variant_loc = parse_variant(variant)
  #print(variant_loc)
  for i in range(len(variant_loc["alts"])):

    # Get the metrics
    var_chrom = variant_loc["chrom"]
    var_start = variant_loc["start"] 
    var_end = variant_loc["end"][i]
    var_type = variant_loc["var_types"][i]
    var_ref = variant_loc["ref"]
    var_alt = variant_loc["alts"][i]
      
    for pileupcolumn in bam.pileup(var_chrom, int(var_start), int(var_end),
                                       min_base_quality = int(min_base_qual),
                                       min_mapping_quality = int(min_mapping_qual),
                                       truncate = True, ignore_overlaps = False):
        
      if pileupcolumn.pos == int(var_start):
        #print(variant_loc)
        summarised_pileup = summarise_pileup(pileupcolumn, var_type, var_alt, var_ref)
        variants_metrics.append(summarised_pileup)
      
  if len(variants_metrics) == 0:
    zero_set = {"ALT": 0 ,"REF": 0, "N": 0, "A": 0, "T": 0, "G": 0, "C": 0, "INS": 0, "DEL": 0, "MNP": 0, "VAR": 0}
    variants_metrics = [{"alt_counts": 0,
                          "ref_counts": 0,
                          "total_counts": 0,
                          "vaf": 0, 
                          "base_counts": zero_set, 
                          "nontrivial_families": zero_set, 
                          "reverse_reads": zero_set, 
                          "umi_discordant": zero_set, 
                          "max_familysize": zero_set, 
                          "mean_bq": zero_set, 
                          "mean_tlen": zero_set}]

  return(variants_metrics)


def update_header(new_header):
  
  new_header.add_meta('FORMAT', items=[('ID',"AD"), ('Number',"A"), ('Type','Integer'), ('Description','Allele depth computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"RD"), ('Number', 1), ('Type','Integer'), ('Description','Reference allele depth computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"DP"), ('Number', 1), ('Type','Integer'), ('Description','Total depth computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"VAF"), ('Number', "A"), ('Type','Float'), ('Description','VAf computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BC_ALT"), ('Number', "A"), ('Type','Integer'), ('Description','Count of ALT bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BC_REF"), ('Number', "A"), ('Type','Integer'), ('Description','Count of REF bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BC_N"), ('Number', "A"), ('Type','Integer'), ('Description','Count of N bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BC_A"), ('Number', "A"), ('Type','Integer'), ('Description','Count of A bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BC_T"), ('Number', "A"), ('Type','Integer'), ('Description','Count of T bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BC_G"), ('Number', "A"), ('Type','Integer'), ('Description','Count of G bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BC_C"), ('Number', "A"), ('Type','Integer'), ('Description','Count of C bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BC_INS"), ('Number', "A"), ('Type','Integer'), ('Description','Count of insertions at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BC_DEL"), ('Number', "A"), ('Type','Integer'), ('Description','Count of deletions bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BC_MNP"), ('Number', "A"), ('Type','Integer'), ('Description','Count of MNPs bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BC_VAR"), ('Number', "A"), ('Type','Integer'), ('Description','Count of other events at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"NTF_ALT"), ('Number', "A"), ('Type','Integer'), ('Description','Count of non-trivial families for ALT bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"NTF_REF"), ('Number', "A"), ('Type','Integer'), ('Description','Count of non-trivial families for REF bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"NTF_N"), ('Number', "A"), ('Type','Integer'), ('Description','Count of non-trivial families for N bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"NTF_A"), ('Number', "A"), ('Type','Integer'), ('Description','Count of non-trivial families for A bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"NTF_T"), ('Number', "A"), ('Type','Integer'), ('Description','Count of non-trivial families for T bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"NTF_G"), ('Number', "A"), ('Type','Integer'), ('Description','Count of non-trivial families for G bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"NTF_C"), ('Number', "A"), ('Type','Integer'), ('Description','Count of non-trivial families for C bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"NTF_INS"), ('Number', "A"), ('Type','Integer'), ('Description','Count of inon-trivial families for nsertions at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"NTF_DEL"), ('Number', "A"), ('Type','Integer'), ('Description','Count of non-trivial families for deletions bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"NTF_MNP"), ('Number', "A"), ('Type','Integer'), ('Description','Count of non-trivial families for MNPs bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"NTF_VAR"), ('Number', "A"), ('Type','Integer'), ('Description','Count of non-trivial families for other events at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"REV_ALT"), ('Number', "A"), ('Type','Integer'), ('Description','Count of ALT bases in reverse reads at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"REV_REF"), ('Number', "A"), ('Type','Integer'), ('Description','Count of REF bases in reverse reads at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"REV_N"), ('Number', "A"), ('Type','Integer'), ('Description','Count of N bases in reverse reads at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"REV_A"), ('Number', "A"), ('Type','Integer'), ('Description','Count of A bases in reverse reads at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"REV_T"), ('Number', "A"), ('Type','Integer'), ('Description','Count of T bases in reverse reads at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"REV_G"), ('Number', "A"), ('Type','Integer'), ('Description','Count of G bases in reverse reads at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"REV_C"), ('Number', "A"), ('Type','Integer'), ('Description','Count of C bases in reverse reads at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"REV_INS"), ('Number', "A"), ('Type','Integer'), ('Description','Count of insertions in reverse reads at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"REV_DEL"), ('Number', "A"), ('Type','Integer'), ('Description','Count of deletions in reverse reads bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"REV_MNP"), ('Number', "A"), ('Type','Integer'), ('Description','Count of MNPs bases in reverse reads at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"REV_VAR"), ('Number', "A"), ('Type','Integer'), ('Description','Count of other events in reverse reads at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"DIS_ALT"), ('Number', "A"), ('Type','Integer'), ('Description','Count of discordant UMI families for ALT bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"DIS_REF"), ('Number', "A"), ('Type','Integer'), ('Description','Count of discordant UMI families for REF bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"DIS_N"), ('Number', "A"), ('Type','Integer'), ('Description','Count of discordant UMI families for N bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"DIS_A"), ('Number', "A"), ('Type','Integer'), ('Description','Count of discordant UMI families for A bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"DIS_T"), ('Number', "A"), ('Type','Integer'), ('Description','Count of discordant UMI families for T bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"DIS_G"), ('Number', "A"), ('Type','Integer'), ('Description','Count of discordant UMI families for G bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"DIS_C"), ('Number', "A"), ('Type','Integer'), ('Description','Count of discordant UMI families for C bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"DIS_INS"), ('Number', "A"), ('Type','Integer'), ('Description','Count of discordant UMI families for insertions at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"DIS_DEL"), ('Number', "A"), ('Type','Integer'), ('Description','Count of discordant UMI families for deletions at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"DIS_MNP"), ('Number', "A"), ('Type','Integer'), ('Description','Count of discordant UMI families for MNPs bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"DIS_VAR"), ('Number', "A"), ('Type','Integer'), ('Description','Count of discordant UMI families for other events at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"MFS_ALT"), ('Number', "A"), ('Type','Float'), ('Description','Mean family size for ALT bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"MFS_REF"), ('Number', "A"), ('Type','Float'), ('Description','Mean family size for REF bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"MFS_N"), ('Number', "A"), ('Type','Float'), ('Description','Mean family size for N bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"MFS_A"), ('Number', "A"), ('Type','Float'), ('Description','Mean family size for A bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"MFS_T"), ('Number', "A"), ('Type','Float'), ('Description','Mean family size for T bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"MFS_G"), ('Number', "A"), ('Type','Float'), ('Description','Mean family size for G bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"MFS_C"), ('Number', "A"), ('Type','Float'), ('Description','Mean family size for C bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"MFS_INS"), ('Number', "A"), ('Type','Float'), ('Description','Mean family size for insertions at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"MFS_DEL"), ('Number', "A"), ('Type','Float'), ('Description','Mean family size for deletions at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"MFS_MNP"), ('Number', "A"), ('Type','Float'), ('Description','Mean family size for MNPs bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"MFS_VAR"), ('Number', "A"), ('Type','Float'), ('Description','Mean family size for other events at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BQ_ALT"), ('Number', "A"), ('Type','Float'), ('Description','Mean base quality for ALT bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BQ_REF"), ('Number', "A"), ('Type','Float'), ('Description','Mean base quality for REF bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BQ_N"), ('Number', "A"), ('Type','Float'), ('Description','Mean base quality for N bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BQ_A"), ('Number', "A"), ('Type','Float'), ('Description','Mean base quality for A bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BQ_T"), ('Number', "A"), ('Type','Float'), ('Description','Mean base quality for T bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BQ_G"), ('Number', "A"), ('Type','Float'), ('Description','Mean base quality for G bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BQ_C"), ('Number', "A"), ('Type','Float'), ('Description','Mean base quality for C bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BQ_INS"), ('Number', "A"), ('Type','Float'), ('Description','Mean base quality for insertions at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BQ_DEL"), ('Number', "A"), ('Type','Float'), ('Description','Mean base quality for deletions at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BQ_MNP"), ('Number', "A"), ('Type','Float'), ('Description','Mean base quality for MNPs bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"BQ_VAR"), ('Number', "A"), ('Type','Float'), ('Description','Mean base quality for other events at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"TL_ALT"), ('Number', "A"), ('Type','Float'), ('Description','Mean fragment length for ALT bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"TL_REF"), ('Number', "A"), ('Type','Float'), ('Description','Mean fragment length for REF bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"TL_N"), ('Number', "A"), ('Type','Float'), ('Description','Mean fragment length for N bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"TL_A"), ('Number', "A"), ('Type','Float'), ('Description','Mean fragment length for A bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"TL_T"), ('Number', "A"), ('Type','Float'), ('Description','Mean fragment length for T bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"TL_G"), ('Number', "A"), ('Type','Float'), ('Description','Mean fragment length for G bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"TL_C"), ('Number', "A"), ('Type','Float'), ('Description','Mean fragment length for C bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"TL_INS"), ('Number', "A"), ('Type','Float'), ('Description','Mean fragment length for insertions at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"TL_DEL"), ('Number', "A"), ('Type','Float'), ('Description','Mean fragment length for deletions at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"TL_MNP"), ('Number', "A"), ('Type','Float'), ('Description','Mean fragment length for MNPs bases at this locus computed by AULE low vaf')])
  new_header.add_meta('FORMAT', items=[('ID',"TL_VAR"), ('Number', "A"), ('Type','Float'), ('Description','Mean fragment length for other events at this locus computed by AULE low vaf')])
  
  return(new_header)

# For making editable entries
def copy_entry(entry, header):
  ret = header.new_record(contig=entry.chrom, start=entry.start, stop=entry.stop,
                          alleles=entry.alleles, id=entry.id, qual=entry.qual)

  return(ret)

# Add formats 
def add_formats(variant, sample_name, variant_metrics):
  
  variant.samples[sample_name]["AD_lv"] = [alt['alt_counts'] for alt in variant_metrics]
  variant.samples[sample_name]["RD_lv"] = [alt['ref_counts'] for alt in variant_metrics]
  variant.samples[sample_name]["DP_lv"] = [alt['total_counts'] for alt in variant_metrics]
  variant.samples[sample_name]["VAF_lv"] = [alt['vaf'] for alt in variant_metrics]
  variant.samples[sample_name]["BC_N"] = [alt['base_counts']["N"] for alt in variant_metrics]
  variant.samples[sample_name]["BC_A"] = [alt['base_counts']["A"] for alt in variant_metrics]
  variant.samples[sample_name]["BC_T"] = [alt['base_counts']["T"] for alt in variant_metrics]
  variant.samples[sample_name]["BC_G"] = [alt['base_counts']["G"] for alt in variant_metrics]
  variant.samples[sample_name]["BC_C"] = [alt['base_counts']["C"] for alt in variant_metrics]
  variant.samples[sample_name]["BC_INS"] = [alt['base_counts']["INS"] for alt in variant_metrics]
  variant.samples[sample_name]["BC_DEL"] = [alt['base_counts']["DEL"] for alt in variant_metrics]
  variant.samples[sample_name]["BC_MNP"] = [alt['base_counts']["MNP"] for alt in variant_metrics]
  variant.samples[sample_name]["BC_VAR"] = [alt['base_counts']["VAR"] for alt in variant_metrics]
  variant.samples[sample_name]["NTF_N"] = [alt['nontrivial_families']["N"] for alt in variant_metrics]
  variant.samples[sample_name]["NTF_A"] = [alt['nontrivial_families']["A"] for alt in variant_metrics]
  variant.samples[sample_name]["NTF_T"] = [alt['nontrivial_families']["T"] for alt in variant_metrics]
  variant.samples[sample_name]["NTF_G"] = [alt['nontrivial_families']["G"] for alt in variant_metrics]
  variant.samples[sample_name]["NTF_C"] = [alt['nontrivial_families']["C"] for alt in variant_metrics]
  variant.samples[sample_name]["NTF_INS"] = [alt['nontrivial_families']["INS"] for alt in variant_metrics]
  variant.samples[sample_name]["NTF_DEL"] = [alt['nontrivial_families']["DEL"] for alt in variant_metrics]
  variant.samples[sample_name]["NTF_MNP"] = [alt['nontrivial_families']["MNP"] for alt in variant_metrics]
  variant.samples[sample_name]["NTF_VAR"] = [alt['nontrivial_families']["VAR"] for alt in variant_metrics]
  variant.samples[sample_name]["REV_N"] = [alt['reverse_reads']["N"] for alt in variant_metrics]
  variant.samples[sample_name]["REV_A"] = [alt['reverse_reads']["A"] for alt in variant_metrics]
  variant.samples[sample_name]["REV_T"] = [alt['reverse_reads']["T"] for alt in variant_metrics]
  variant.samples[sample_name]["REV_G"] = [alt['reverse_reads']["G"] for alt in variant_metrics]
  variant.samples[sample_name]["REV_C"] = [alt['reverse_reads']["C"] for alt in variant_metrics]
  variant.samples[sample_name]["REV_INS"] = [alt['reverse_reads']["INS"] for alt in variant_metrics]
  variant.samples[sample_name]["REV_DEL"] = [alt['reverse_reads']["DEL"] for alt in variant_metrics]
  variant.samples[sample_name]["REV_MNP"] = [alt['reverse_reads']["MNP"] for alt in variant_metrics]
  variant.samples[sample_name]["REV_VAR"] = [alt['reverse_reads']["VAR"] for alt in variant_metrics]
  variant.samples[sample_name]["DIS_N"] = [alt['umi_discordant']["N"] for alt in variant_metrics]
  variant.samples[sample_name]["DIS_A"] = [alt['umi_discordant']["A"] for alt in variant_metrics]
  variant.samples[sample_name]["DIS_T"] = [alt['umi_discordant']["T"] for alt in variant_metrics]
  variant.samples[sample_name]["DIS_G"] = [alt['umi_discordant']["G"] for alt in variant_metrics]
  variant.samples[sample_name]["DIS_C"] = [alt['umi_discordant']["C"] for alt in variant_metrics]
  variant.samples[sample_name]["DIS_INS"] = [alt['umi_discordant']["INS"] for alt in variant_metrics]
  variant.samples[sample_name]["DIS_DEL"] = [alt['umi_discordant']["DEL"] for alt in variant_metrics]
  variant.samples[sample_name]["DIS_MNP"] = [alt['umi_discordant']["MNP"] for alt in variant_metrics]
  variant.samples[sample_name]["DIS_VAR"] = [alt['umi_discordant']["VAR"] for alt in variant_metrics]
  variant.samples[sample_name]["MFS_N"] = [alt['max_familysize']["N"] for alt in variant_metrics]
  variant.samples[sample_name]["MFS_A"] = [alt['max_familysize']["A"] for alt in variant_metrics]
  variant.samples[sample_name]["MFS_T"] = [alt['max_familysize']["T"] for alt in variant_metrics]
  variant.samples[sample_name]["MFS_G"] = [alt['max_familysize']["G"] for alt in variant_metrics]
  variant.samples[sample_name]["MFS_C"] = [alt['max_familysize']["C"] for alt in variant_metrics]
  variant.samples[sample_name]["MFS_INS"] = [alt['max_familysize']["INS"] for alt in variant_metrics]
  variant.samples[sample_name]["MFS_DEL"] = [alt['max_familysize']["DEL"] for alt in variant_metrics]
  variant.samples[sample_name]["MFS_MNP"] = [alt['max_familysize']["MNP"] for alt in variant_metrics]
  variant.samples[sample_name]["MFS_VAR"] = [alt['max_familysize']["VAR"] for alt in variant_metrics]
  variant.samples[sample_name]["BQ_N"] = [alt['mean_bq']["N"] for alt in variant_metrics]
  variant.samples[sample_name]["BQ_A"] = [alt['mean_bq']["A"] for alt in variant_metrics]
  variant.samples[sample_name]["BQ_T"] = [alt['mean_bq']["T"] for alt in variant_metrics]
  variant.samples[sample_name]["BQ_G"] = [alt['mean_bq']["G"] for alt in variant_metrics]
  variant.samples[sample_name]["BQ_C"] = [alt['mean_bq']["C"] for alt in variant_metrics]
  variant.samples[sample_name]["BQ_INS"] = [alt['mean_bq']["INS"] for alt in variant_metrics]
  variant.samples[sample_name]["BQ_DEL"] = [alt['mean_bq']["DEL"] for alt in variant_metrics]
  variant.samples[sample_name]["BQ_MNP"] = [alt['mean_bq']["MNP"] for alt in variant_metrics]
  variant.samples[sample_name]["BQ_VAR"] = [alt['mean_bq']["VAR"] for alt in variant_metrics]
  variant.samples[sample_name]["TL_N"] = [alt['mean_tlen']["N"] for alt in variant_metrics]
  variant.samples[sample_name]["TL_A"] = [alt['mean_tlen']["A"] for alt in variant_metrics]
  variant.samples[sample_name]["TL_T"] = [alt['mean_tlen']["T"] for alt in variant_metrics]
  variant.samples[sample_name]["TL_G"] = [alt['mean_tlen']["G"] for alt in variant_metrics]
  variant.samples[sample_name]["TL_C"] = [alt['mean_tlen']["C"] for alt in variant_metrics]
  variant.samples[sample_name]["TL_INS"] = [alt['mean_tlen']["INS"] for alt in variant_metrics]
  variant.samples[sample_name]["TL_DEL"] = [alt['mean_tlen']["DEL"] for alt in variant_metrics]
  variant.samples[sample_name]["TL_MNP"] = [alt['mean_tlen']["MNP"] for alt in variant_metrics]
  variant.samples[sample_name]["TL_VAR"] = [alt['mean_tlen']["VAR"] for alt in variant_metrics]
  
  return(variant)


def main():

    # Parse arguments 
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf_in', metavar='<vcf_in>', dest="vcf_in", help="Input vcf")
    parser.add_argument('-vcf_out', metavar='<vcf_out>', dest="vcf_out", help="Annotated output VCF")
    parser.add_argument('-tumour_bam', metavar='<tumour_bam>', dest="tumour_bam", help="Input tumour bam file")
    parser.add_argument('-normal_bam', metavar='<normal_bam>', dest="normal_bam", help="Input matched normal bam file")
    parser.add_argument('-threads', metavar='<threads>', dest="threads", help="Number of parallel threads")
    parser.add_argument('-genome_build', metavar='<genome_build>', dest="genome_build", help="Genome build (hg19 or hg38)")
    parser.add_argument('-min_base_qual', metavar='<min_base_qual>', dest="min_base_qual", help="Minimum base quality")
    parser.add_argument('-min_mapping_qual', metavar='<min_mapping_qual>', dest="min_mapping_qual", help="Minimum MAPQ")

    args = parser.parse_args()
    
    # Load the VCF file
    vcf = pysam.VariantFile(args.vcf_in)
    
    # Load BAM files 
    tumour_bam = pysam.AlignmentFile(args.tumour_bam)
    normal_bam = pysam.AlignmentFile(args.normal_bam)
    
    # Update a header for a new VCF file 
    new_header =  vcf.header.copy()
    new_header =  update_header(new_header)
    
    # Add new samples to the vcf file 
    new_header.add_sample("TUMOUR_aule")
    new_header.add_sample("NORMAL_aule")
    
    # Initiate output vcf 
    with pysam.VariantFile(args.vcf_out, "w", header=new_header) as out_vcf:
   # out_vcf = pysam.VariantFile(args.vcf_out, 'w', header=new_header)
    
      for variant in vcf.fetch():
        
        updated_variant = copy_entry(variant, new_header)
        
        # Extract tumour sample metrics
        variant_metrics_tumour = get_variant_metrics(updated_variant, tumour_bam, args.min_base_qual, args.min_mapping_qual, new_header)
        #updated_variant = add_formats(updated_variant, "TUMOUR_aule", variant_metrics_tumour)
        #print("TUMOUR:", variant_metrics_tumour)
        # Extract normal sample metrics
        variant_metrics_normal = get_variant_metrics(updated_variant, normal_bam, args.min_base_qual, args.min_mapping_qual, new_header)
        #updated_variant = add_formats(updated_variant, "NORMAL_aule", variant_metrics_normal) # It doesn't work insode the function
        #print("NORMAL:", variant_metrics_normal)
        
        # Assign values - tumour sample
        if len(variant_metrics_tumour) > 0:
          updated_variant.samples["TUMOUR_aule"]["AD"] = [alt['alt_counts'] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["RD"] = [alt['ref_counts'] for alt in variant_metrics_tumour][0]
          updated_variant.samples["TUMOUR_aule"]["DP"] = [alt['total_counts'] for alt in variant_metrics_tumour][0]
          updated_variant.samples["TUMOUR_aule"]["VAF"] = [alt['vaf'] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BC_ALT"] = [alt['base_counts']["ALT"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BC_REF"] = [alt['base_counts']["REF"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["NTF_ALT"] = [alt['nontrivial_families']["ALT"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["NTF_REF"] = [alt['nontrivial_families']["REF"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["REV_ALT"] = [alt['reverse_reads']["ALT"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["REV_REF"] = [alt['reverse_reads']["REF"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["DIS_ALT"] = [alt['umi_discordant']["ALT"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["DIS_REF"] = [alt['umi_discordant']["REF"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["MFS_ALT"] = [alt['max_familysize']["ALT"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["MFS_REF"] = [alt['max_familysize']["REF"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BQ_ALT"] = [alt['mean_bq']["ALT"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BQ_REF"] = [alt['mean_bq']["REF"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["TL_ALT"] = [alt['mean_tlen']["ALT"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["TL_REF"] = [alt['mean_tlen']["REF"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BC_N"] = [alt['base_counts']["N"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BC_A"] = [alt['base_counts']["A"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BC_T"] = [alt['base_counts']["T"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BC_G"] = [alt['base_counts']["G"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BC_C"] = [alt['base_counts']["C"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BC_INS"] = [alt['base_counts']["INS"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BC_DEL"] = [alt['base_counts']["DEL"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BC_MNP"] = [alt['base_counts']["MNP"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BC_VAR"] = [alt['base_counts']["VAR"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["NTF_N"] = [alt['nontrivial_families']["N"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["NTF_A"] = [alt['nontrivial_families']["A"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["NTF_T"] = [alt['nontrivial_families']["T"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["NTF_G"] = [alt['nontrivial_families']["G"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["NTF_C"] = [alt['nontrivial_families']["C"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["NTF_INS"] = [alt['nontrivial_families']["INS"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["NTF_DEL"] = [alt['nontrivial_families']["DEL"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["NTF_MNP"] = [alt['nontrivial_families']["MNP"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["NTF_VAR"] = [alt['nontrivial_families']["VAR"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["REV_N"] = [alt['reverse_reads']["N"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["REV_A"] = [alt['reverse_reads']["A"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["REV_T"] = [alt['reverse_reads']["T"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["REV_G"] = [alt['reverse_reads']["G"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["REV_C"] = [alt['reverse_reads']["C"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["REV_INS"] = [alt['reverse_reads']["INS"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["REV_DEL"] = [alt['reverse_reads']["DEL"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["REV_MNP"] = [alt['reverse_reads']["MNP"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["REV_VAR"] = [alt['reverse_reads']["VAR"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["DIS_N"] = [alt['umi_discordant']["N"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["DIS_A"] = [alt['umi_discordant']["A"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["DIS_T"] = [alt['umi_discordant']["T"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["DIS_G"] = [alt['umi_discordant']["G"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["DIS_C"] = [alt['umi_discordant']["C"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["DIS_INS"] = [alt['umi_discordant']["INS"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["DIS_DEL"] = [alt['umi_discordant']["DEL"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["DIS_MNP"] = [alt['umi_discordant']["MNP"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["DIS_VAR"] = [alt['umi_discordant']["VAR"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["MFS_N"] = [alt['max_familysize']["N"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["MFS_A"] = [alt['max_familysize']["A"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["MFS_T"] = [alt['max_familysize']["T"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["MFS_G"] = [alt['max_familysize']["G"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["MFS_C"] = [alt['max_familysize']["C"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["MFS_INS"] = [alt['max_familysize']["INS"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["MFS_DEL"] = [alt['max_familysize']["DEL"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["MFS_MNP"] = [alt['max_familysize']["MNP"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["MFS_VAR"] = [alt['max_familysize']["VAR"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BQ_N"] = [alt['mean_bq']["N"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BQ_A"] = [alt['mean_bq']["A"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BQ_T"] = [alt['mean_bq']["T"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BQ_G"] = [alt['mean_bq']["G"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BQ_C"] = [alt['mean_bq']["C"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BQ_INS"] = [alt['mean_bq']["INS"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BQ_DEL"] = [alt['mean_bq']["DEL"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BQ_MNP"] = [alt['mean_bq']["MNP"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["BQ_VAR"] = [alt['mean_bq']["VAR"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["TL_N"] = [alt['mean_tlen']["N"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["TL_A"] = [alt['mean_tlen']["A"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["TL_T"] = [alt['mean_tlen']["T"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["TL_G"] = [alt['mean_tlen']["G"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["TL_C"] = [alt['mean_tlen']["C"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["TL_INS"] = [alt['mean_tlen']["INS"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["TL_DEL"] = [alt['mean_tlen']["DEL"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["TL_MNP"] = [alt['mean_tlen']["MNP"] for alt in variant_metrics_tumour]
          updated_variant.samples["TUMOUR_aule"]["TL_VAR"] = [alt['mean_tlen']["VAR"] for alt in variant_metrics_tumour]
          
          # Assign values - normal sample
          #print("NORM", variant_metrics_normal)
          updated_variant.samples["NORMAL_aule"]["AD"] = [alt['alt_counts'] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["RD"] = [alt['ref_counts'] for alt in variant_metrics_normal][0]
          updated_variant.samples["NORMAL_aule"]["DP"] = [alt['total_counts'] for alt in variant_metrics_normal][0]
          updated_variant.samples["NORMAL_aule"]["VAF"] = [alt['vaf'] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BC_ALT"] = [alt['base_counts']["ALT"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BC_REF"] = [alt['base_counts']["REF"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["NTF_ALT"] = [alt['nontrivial_families']["ALT"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["NTF_REF"] = [alt['nontrivial_families']["REF"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["REV_ALT"] = [alt['reverse_reads']["ALT"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["REV_REF"] = [alt['reverse_reads']["REF"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["DIS_ALT"] = [alt['umi_discordant']["ALT"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["DIS_REF"] = [alt['umi_discordant']["REF"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["MFS_ALT"] = [alt['max_familysize']["ALT"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["MFS_REF"] = [alt['max_familysize']["REF"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BQ_ALT"] = [alt['mean_bq']["ALT"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BQ_REF"] = [alt['mean_bq']["REF"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["TL_ALT"] = [alt['mean_tlen']["ALT"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["TL_REF"] = [alt['mean_tlen']["REF"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BC_N"] = [alt['base_counts']["N"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BC_A"] = [alt['base_counts']["A"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BC_T"] = [alt['base_counts']["T"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BC_G"] = [alt['base_counts']["G"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BC_C"] = [alt['base_counts']["C"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BC_INS"] = [alt['base_counts']["INS"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BC_DEL"] = [alt['base_counts']["DEL"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BC_MNP"] = [alt['base_counts']["MNP"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BC_VAR"] = [alt['base_counts']["VAR"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["NTF_N"] = [alt['nontrivial_families']["N"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["NTF_A"] = [alt['nontrivial_families']["A"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["NTF_T"] = [alt['nontrivial_families']["T"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["NTF_G"] = [alt['nontrivial_families']["G"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["NTF_C"] = [alt['nontrivial_families']["C"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["NTF_INS"] = [alt['nontrivial_families']["INS"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["NTF_DEL"] = [alt['nontrivial_families']["DEL"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["NTF_MNP"] = [alt['nontrivial_families']["MNP"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["NTF_VAR"] = [alt['nontrivial_families']["VAR"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["REV_N"] = [alt['reverse_reads']["N"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["REV_A"] = [alt['reverse_reads']["A"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["REV_T"] = [alt['reverse_reads']["T"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["REV_G"] = [alt['reverse_reads']["G"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["REV_C"] = [alt['reverse_reads']["C"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["REV_INS"] = [alt['reverse_reads']["INS"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["REV_DEL"] = [alt['reverse_reads']["DEL"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["REV_MNP"] = [alt['reverse_reads']["MNP"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["REV_VAR"] = [alt['reverse_reads']["VAR"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["DIS_N"] = [alt['umi_discordant']["N"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["DIS_A"] = [alt['umi_discordant']["A"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["DIS_T"] = [alt['umi_discordant']["T"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["DIS_G"] = [alt['umi_discordant']["G"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["DIS_C"] = [alt['umi_discordant']["C"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["DIS_INS"] = [alt['umi_discordant']["INS"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["DIS_DEL"] = [alt['umi_discordant']["DEL"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["DIS_MNP"] = [alt['umi_discordant']["MNP"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["DIS_VAR"] = [alt['umi_discordant']["VAR"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["MFS_N"] = [alt['max_familysize']["N"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["MFS_A"] = [alt['max_familysize']["A"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["MFS_T"] = [alt['max_familysize']["T"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["MFS_G"] = [alt['max_familysize']["G"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["MFS_C"] = [alt['max_familysize']["C"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["MFS_INS"] = [alt['max_familysize']["INS"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["MFS_DEL"] = [alt['max_familysize']["DEL"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["MFS_MNP"] = [alt['max_familysize']["MNP"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["MFS_VAR"] = [alt['max_familysize']["VAR"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BQ_N"] = [alt['mean_bq']["N"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BQ_A"] = [alt['mean_bq']["A"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BQ_T"] = [alt['mean_bq']["T"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BQ_G"] = [alt['mean_bq']["G"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BQ_C"] = [alt['mean_bq']["C"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BQ_INS"] = [alt['mean_bq']["INS"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BQ_DEL"] = [alt['mean_bq']["DEL"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BQ_MNP"] = [alt['mean_bq']["MNP"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["BQ_VAR"] = [alt['mean_bq']["VAR"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["TL_N"] = [alt['mean_tlen']["N"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["TL_A"] = [alt['mean_tlen']["A"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["TL_T"] = [alt['mean_tlen']["T"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["TL_G"] = [alt['mean_tlen']["G"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["TL_C"] = [alt['mean_tlen']["C"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["TL_INS"] = [alt['mean_tlen']["INS"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["TL_DEL"] = [alt['mean_tlen']["DEL"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["TL_MNP"] = [alt['mean_tlen']["MNP"] for alt in variant_metrics_normal]
          updated_variant.samples["NORMAL_aule"]["TL_VAR"] = [alt['mean_tlen']["VAR"] for alt in variant_metrics_normal]
          
        out_vcf.write(updated_variant)
      
    out_vcf.close()
   
    print("Done!")

if __name__ == "__main__":
	main()
