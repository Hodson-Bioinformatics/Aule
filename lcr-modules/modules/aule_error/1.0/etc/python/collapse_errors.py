import argparse
import os
import numpy as np
import pandas as pd
import pysam
import time

def main():
  global args
  
  parser = argparse.ArgumentParser()
  parser.add_argument('-in_bam', metavar='<in_bam>', dest="in_bam", help="input bam file")
  # Input bam is output of CallMolecularConsensusReads, re-aligned
  parser.add_argument('-in_coords', metavar='<in_coords>', dest="in_coords", help="input coords file")
  parser.add_argument('-out_errors', metavar='<out_errors>', dest="out_errors", help="output errors csv file")
  parser.add_argument('-out_stats', metavar='<out_stats>', dest="out_stats", help="output stats csv file")
  
  args = parser.parse_args()
  
  collapse_errors(args.in_bam, args.in_coords, args.out_errors, args.out_stats)

def collapse_errors(in_bam, in_coords, out_errors, out_stats):
  coords = pd.read_csv(in_coords, encoding = "latin-1")
  reverse_lookup = {x:i for i, x in enumerate(list(coords["chr_pos"].values))}
  chr_key = {i:s for i,s in enumerate([f"chr{x}" for x in list(range(1, 23))+["X","Y"]], 0)}
  
  bam = pysam.Samfile(in_bam, "rb")
  
  stats = pd.DataFrame({"reads": [0], "families": [0], "families_nontrivial": [0], "families_size1": [0], "families_size2": [0], 
  "families_size3": [0], "total_fragment_size": [0], "fragment_size_on_target": [0], "fragments_on_target": [0], "time_elapsed": [0]})
  
  start_time = time.time()
  for read in bam.fetch(until_eof = True):
    if read.is_read1:
      read1 = read
    else:
      read2 = read
      
      if read1.rname <= 23 and read1.rname >= 0:
        # coords, total_fg_size, fg_size_aln_with_panel, n_fgs_on_target = collapse_family_info(coords, read1, read2, total_fg_size, fg_size_aln_with_panel, n_fgs_on_target, reverse_lookup)
        coords, stats = collapse_family_info(coords, stats, read1, read2, reverse_lookup, chr_key)
        # stats.loc[0, "time_elapsed"] = time.time() - start_time
        stats.iloc[0, 9] = time.time() - start_time
      
      ### temporary testing settings
      # print_stats(stats, every_n = 1000)
      # if stats["families"][0] >= 10000:
      #   break
      print_stats(stats, every_n = 1000000, file = out_stats)
      
  ### add a bit to tidy both dfs before saving them
  # coords["error_rate"], coords["error_rate_expanded"] = coords.apply(add_statistics_per_position, axis = 1)
  
  coords.to_csv(out_errors, index = False)
  stats.to_csv(out_stats, index = False)
  
  print("----------------------------------------------------------------------------------------------------")
  print("Totals:")
  print("     -", stats["reads"][0], "from", stats["families"][0], "UMI families, of which")
  print("     -", stats["families_nontrivial"][0], "were of size 2 or greater (" + str(stats["families_nontrivial"][0] * 100 // stats["families"][0]) + "%)")
  print("     -", stats["fragments_on_target"][0], "fragments on target (" + str(round(stats["fragments_on_target"][0] * 100 / stats["families"][0])) + "%)")
  print("     -", str(round(stats["fragment_size_on_target"][0] * 100 / stats["total_fragment_size"][0])) + "%", "of total fragment length aligns with panel")
  print("     -", "time elapsed:", pretty_time_delta(stats["time_elapsed"][0]))
  print("Information was saved in " + out_errors)
  print("----------------------------------------------------------------------------------------------------")

# replaced by separate script, outputting the coords df in csv
def bed_to_coords(bed):
  # sort
  chr_key={s:i for i,s in enumerate([f"chr{x}" for x in list(range(1,23))+["X","Y"]],1)}
  bed["chrid"] = bed["chr"].map(chr_key)
  bed = bed.sort_values(["chrid", "start", "end"], ascending = True).reset_index(drop=True)
  
  # expand_regions
  for index, row in bed.iterrows():
    temp = pd.DataFrame({"name": row["name"], "chr": row["chr"], "pos": range(row["start"], row["end"] + 1)})
    if index == 0:
      coords = temp
    else:
      coords = pd.concat([coords, temp])
  
  coords.reset_index(inplace = True)
  coords["index"] = coords.index
  
  coords["chr_pos"] = [chr + "_" + str(pos) for chr, pos in zip(coords["chr"], coords["pos"])]
  
  # add cols
  coords["reads"] = 0
  coords["families"] = 0
  coords["reads_nontrivial"] = 0
  coords["families_nontrivial"] = 0
  coords["families_size1"] = 0
  coords["families_size2"] = 0
  coords["families_size3"] = 0
  coords["r12_overlap_depth"] = 0
  coords["total_base_quality"] = 0
  coords["total_mapping_quality"] = 0
  coords["total_distance_from_5p"] = 0
  coords["total_distance_from_3p"] = 0
  coords["total_fragment_length"] = 0
  coords["umi_errors"] = 0
  coords["r12_discordant"] = 0
  coords["aligned_base_A"] = 0
  coords["aligned_base_T"] = 0
  coords["aligned_base_C"] = 0
  coords["aligned_base_G"] = 0
  coords["aligned_base_nontrivial_A"] = 0
  coords["aligned_base_nontrivial_T"] = 0
  coords["aligned_base_nontrivial_C"] = 0
  coords["aligned_base_nontrivial_G"] = 0
  
  return coords

# @profile
def collapse_family_info(coords, stats, read1, read2, reverse_lookup, chr_key):
  # if reverse_lookup is None:
  #   reverse_lookup = {x:i for i, x in enumerate(list(coords["chr_pos"].values))}
  # 
  # if chr_key is None:
  #   chr_key = {i:s for i,s in enumerate([f"chr{x}" for x in list(range(1,23))+["X","Y"]],0)}
  
  chrom = chr_key[read1.rname]
  
  # get the indices of the read that map to ref coords using read.get_aligned_pairs()
  r1_idx = [x[0] for x in read1.get_aligned_pairs(matches_only = True)]
  
  # get cd and ce tags corresponding to ref coords
  r1_cd = np.array(read1.get_tag("cd"))[r1_idx]
  r1_ce = np.array(read1.get_tag("ce"))[r1_idx]
  
  # get corresponding reference positions
  r1_ref = read1.get_reference_positions()
  
  # same for read2 processed separately
  r2_idx = [x[0] for x in read2.get_aligned_pairs(matches_only = True)]
  r2_cd  = np.array(read2.get_tag("cd"))[r2_idx]
  r2_ce  = np.array(read2.get_tag("ce"))[r2_idx]
  r2_ref = read2.get_reference_positions()
  
  # also process fragment as one
  fg_ref = sorted(list(set(r1_ref).union(set(r2_ref)))) # fragment ref positions
    # note these use r1 information for the overlap
  fg_cd = [r1_cd[r1_ref.index(x)] if x in r1_ref else r2_cd[r2_ref.index(x)] for x in fg_ref]
  fg_ce = [r1_ce[r1_ref.index(x)] if x in r1_ref else r2_ce[r2_ref.index(x)] for x in fg_ref]
    # get base qualities of fragment
  fg_qual = [read1.qual[r1_ref.index(x)] if x in r1_ref else read2.qual[r2_ref.index(x)] for x in fg_ref]
  fg_seq = [read1.seq[r1_ref.index(x)] if x in r1_ref else read2.seq[r2_ref.index(x)] for x in fg_ref]
  
  # also determine the discordance between r1 and r2 as part of this
  # (by checking only the overlap bases and seeing if they match)
  overlap_ref = sorted(list(set(r1_ref).intersection(set(r2_ref))))
  r1_overlap_idx = [x[0] for x in read1.get_aligned_pairs(matches_only = True) if x[1] in overlap_ref]
  r1_overlap_seq = np.array([*read1.seq])[r1_overlap_idx]
  r2_overlap_idx = [x[0] for x in read2.get_aligned_pairs(matches_only = True) if x[1] in overlap_ref]
  r2_overlap_seq = np.array([*read2.seq])[r2_overlap_idx]
  overlap_discordance = [0 if r1_overlap_seq[i] == r2_overlap_seq[i] else 1 for i in range(len(overlap_ref))]
  
  # find rows in coords table and add new information
    # first need to subset all of the above to coords that are part of the panel
  def find_indices_in_bed_coords(chrom, ref, reverse_lookup):
    query_chr_pos = [chrom + "_" + str(x + 1) for x in ref] # +1 because ref is 0-indexed and bed is 1-indexed
    return [i for i, x in enumerate(query_chr_pos) if reverse_lookup.get(x, None) is not None]
  
  fg_sub_idx  = find_indices_in_bed_coords(chrom, fg_ref, reverse_lookup)
  fg_ref_sub  = np.array(fg_ref)[fg_sub_idx]  # subset
  fg_cd_sub   = np.array(fg_cd)[fg_sub_idx]   # subset
  fg_ce_sub   = np.array(fg_ce)[fg_sub_idx]   # subset
  fg_qual_sub = np.array(fg_qual)[fg_sub_idx] # subset
  fg_seq_sub  = np.array(fg_seq)[fg_sub_idx]  # subset
  
  overlap_sub_idx         = find_indices_in_bed_coords(chrom, overlap_ref, reverse_lookup)
  overlap_ref_sub         = np.array(overlap_ref)[overlap_sub_idx]         # subset
  overlap_discordance_sub = np.array(overlap_discordance)[overlap_sub_idx] # subset
  
    # for fg and overlap separately, create index and update values
    # ensure the rows of coords are returned in the same order as fg/overlap coords
  coords_fg_idx = [reverse_lookup.get(chrom + "_" + str(x + 1)) for x in fg_ref_sub]
  coords_overlap_idx = [reverse_lookup.get(chrom + "_" + str(x + 1)) for x in overlap_ref_sub]
  
  stats.iloc[0, 0]     += read1.get_tag("cD") # "reads"
  stats.iloc[0, 1]     += 1 # "families"
  if read1.get_tag("cD") == 1:
    stats.iloc[0, 3]   += 1 # "families_size1"
  else:
    stats.iloc[0, 2]   += 1 # "families_nontrivial"
    if read1.get_tag("cD") == 2:
      stats.iloc[0, 4] += 1 # "families_size2"
    elif read1.get_tag("cD") == 3:
      stats.iloc[0, 5] += 1 # "families_size3"
  stats.iloc[0, 6]     += len(fg_ref) # "total_fragment_size"
  stats.iloc[0, 7]     += len(fg_ref_sub) # "fragment_size_on_target"
  stats.iloc[0, 8]     += int(len(fg_ref_sub) > 0) # "fragments_on_target"
  
  coords.iloc[coords_fg_idx, 10]      += fg_cd_sub # "reads"
  coords.iloc[coords_fg_idx, 11]      += 1 # "families"
  coords.iloc[coords_overlap_idx, 17] += 1 # "r12_overlap_depth"
  coords.iloc[coords_fg_idx, 18]      += [ord(x) - 33 for x in fg_qual_sub] # convert ASCII to integers # "total_base_quality"
  coords.iloc[coords_fg_idx, 19]      += read1.mapq # "total_mapping_quality"
  coords.iloc[coords_fg_idx, 20]      += [x + 1 for x in fg_sub_idx] # "total_distance_from_5p"
  coords.iloc[coords_fg_idx, 21]      += [max(fg_sub_idx) + 1 - x for x in fg_sub_idx] # "total_distance_from_3p"
  coords.iloc[coords_fg_idx, 22]      += abs(read1.tlen) # "total_fragment_length"
  coords.iloc[coords_overlap_idx, 24] += overlap_discordance_sub # "r12_discordant"
  coords.iloc[coords_fg_idx, 25]      += [1 if x == "A" else 0 for x in fg_seq_sub] # "aligned_base_A"
  coords.iloc[coords_fg_idx, 26]      += [1 if x == "T" else 0 for x in fg_seq_sub] # "aligned_base_T"
  coords.iloc[coords_fg_idx, 27]      += [1 if x == "C" else 0 for x in fg_seq_sub] # "aligned_base_C"
  coords.iloc[coords_fg_idx, 28]      += [1 if x == "G" else 0 for x in fg_seq_sub] # "aligned_base_G"
  
  if read1.get_tag("cD") == 1:
    coords.iloc[coords_fg_idx, 14]    += 1 # "families_size1"
  elif read1.get_tag("cD") > 1:
    coords.iloc[coords_fg_idx, 12]    += [x if x > 1 else 0 for x in fg_cd_sub] # "reads_nontrivial"
    coords.iloc[coords_fg_idx, 13]    += [1 if x > 1 else 0 for x in fg_cd_sub] # "families_nontrivial"
    coords.iloc[coords_fg_idx, 15]    += [1 if x == 2 else 0 for x in fg_cd_sub] # "families_size2"
    coords.iloc[coords_fg_idx, 16]    += [1 if x == 3 else 0 for x in fg_cd_sub] # "families_size3"
    coords.iloc[coords_fg_idx, 23]    += fg_ce_sub # "umi_errors"
    coords.iloc[coords_fg_idx, 29]    += [1 if (x > 1) & (seq == "A") else 0 for x, seq in zip(fg_cd_sub, fg_seq_sub)] # "aligned_base_nontrivial_A"
    coords.iloc[coords_fg_idx, 30]    += [1 if (x > 1) & (seq == "T") else 0 for x, seq in zip(fg_cd_sub, fg_seq_sub)] # "aligned_base_nontrivial_T"
    coords.iloc[coords_fg_idx, 31]    += [1 if (x > 1) & (seq == "C") else 0 for x, seq in zip(fg_cd_sub, fg_seq_sub)] # "aligned_base_nontrivial_C"
    coords.iloc[coords_fg_idx, 32]    += [1 if (x > 1) & (seq == "G") else 0 for x, seq in zip(fg_cd_sub, fg_seq_sub)] # "aligned_base_nontrivial_G"
  
  return coords, stats

def print_stats(stats, every_n = 1000000, file = None):
  if stats["families"][0] % every_n == 0:
    if file is not None:
      stats.to_csv(file, index = False)
    print(stats["families"][0], "UMI families done,", 
    stats["families_nontrivial"][0], "(" + str(round(stats["families_nontrivial"][0] * 100 / stats["families"][0])) + "%) of size 2 or more;",
    pretty_time_delta(stats["time_elapsed"][0]), "elapsed;", 
    stats["fragments_on_target"][0], "fragments on target",
    "(" + str(round(stats["fragments_on_target"][0] * 100 / stats["families"][0])) + "%)", 
    "-", str(round(stats["fragment_size_on_target"][0] * 100 / stats["total_fragment_size"][0])) + 
    "% of total fragment length aligns with panel")

def pretty_time_delta(seconds):
  sign_string = '-' if seconds < 0 else ''
  seconds = abs(int(seconds))
  years, seconds = divmod(seconds, 31536000)
  days, seconds = divmod(seconds, 86400)
  hours, seconds = divmod(seconds, 3600)
  minutes, seconds = divmod(seconds, 60)
  if years > 0:
    return '%s%dy%dd%dh%dm%ds' % (sign_string, years, days, hours, minutes, seconds)
  elif days > 0:
    return '%s%dd%dh%dm%ds' % (sign_string, days, hours, minutes, seconds)
  elif hours > 0:
    return '%s%dh%dm%ds' % (sign_string, hours, minutes, seconds)
  elif minutes > 0:
    return '%s%dm%ds' % (sign_string, minutes, seconds)
  else:
    return '%s%ds' % (sign_string, seconds)

if __name__ == "__main__":
  main()
