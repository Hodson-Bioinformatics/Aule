#!/usr/bin/env python
# Inputs (via Snakemake):
# 1. one tumour bam (or cram) file or both a tumour and normal bam file
# 2. an index MAF file that will be updated with missing rows. 
# The index MAF should have variants called from analysis of the bam file(s) supplied
# 3. One ore more additional MAFs with variants called from the same patient using other bam files
# Output: a new MAF file that contains all variants from the index MAF plus rows for every variant unique to the additional MAFs
# Important: this script will determine the tumour_ref_count and tumour_alt_count for new MAF rows using the supplied bam
# It uses the pysam pileup engine to calculate the VAF of SNVs, indels based on the bam (i.e. trusting the alignment as-is)

#If you are using tabular data that isn't a MAF, the minimum set of required columns are:
#Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Type, t_ref_count, t_alt_count, t_depth
# e.g. chr2 1234 1235 A T SNP 0 0 0
from pathlib import Path 
import pysam
import os
import pandas as pd
import multiprocessing
import oncopipe as op
from multiprocessing import set_start_method
verbose = True #turn off print statements
min_map_qual = 1 #minimum mapping quality of a read to use it
padding = 100

#restrict processing to only one chromosome, if desired. For parallelization.
#this list should be automatically pulled from the MAFs but I worry that if one or more chromosomes is missing from the index MAF this might lead to issues

#normal_bam = ""


#TODO: Normal bam is not implemented yet! Currently this just extracts read counts for the mutation site in the tumour bam
#augmented_maf = "./00-14595_tumorC.grch37.augmented.maf"
#tumour_bam_file = "data/genome_bams/00-14595_tumorC.grch37.bam"
#normal_bam_file = "data/genome_bams/00-14595_normal.grch37.bam"
#index_maf_file = "results/gambl/slms-3_vcf2maf_current/99-outputs/genome--grch37/00-14595_tumorC--00-14595_normal--matched_slms-3.final.maf"

#add_maf_files = ["results/gambl/slms-3_vcf2maf_current/99-outputs/genome--grch37/00-14595_tumorA--00-14595_normal--matched_slms-3.final.maf",
#                "results/gambl/slms-3_vcf2maf_current/99-outputs/genome--grch37/00-14595_tumorB--00-14595_normal--matched_slms-3.final.maf"]
#remove the above lines and uncomment the next few lines to get this snakemake-ready
nproc = snakemake.threads
genome_build = snakemake.wildcards["genome_build"]

augmented_maf = snakemake.output["augmented_maf"] #your new MAF complete with all variants in index_maf plus the variants unique to all of mafs_to_add
tumour_bam_file = snakemake.input["tumour_bam"] #the bam file you want to check for read support for variants
index_maf_file = snakemake.input["index_maf"] #the MAF you are augmenting (with variants called from tumour_bam)
add_maf_files = snakemake.input["add_maf_files"] #list of paths to the other MAF files you want to fill in. 
normal_bam_file = []

# This checks if patient has a single tumour timepoint. If yes then it will just symlink the slms_3 maf instead of running augment_ssm
if len(add_maf_files) == 0:
    print(f"I'll make a symlink")
    op.absolute_symlink(index_maf_file, augmented_maf)
    exit()
    
tumour_sample_id = snakemake.wildcards["tumour_id"]

#special case: if index_maf_file is among add_maf files, zero all t_ref_count and t_alt_count values
if genome_build in ["hg38","hg19-reddy"]:
    chromosomes = ["chrY","chrX","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]
else:
    chromosomes = ["Y","X","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]



#normal_bam_file = snakemake.input["normal_bam"]
# I hope this will work but I've never tried it with >1 input file so this may need to be tweaked



#try:
#    normal_bam = snakemake.input["normal_bam"]
#except AttributeError:
#    print("proceeding without normal bam file")

# Input: two data frames from MAF files, an index maf and a supplemental MAF
# Output: just the rows from the supplemental MAF that are not present in the index MAF
def get_missing_maf_rows(this_index_maf,other_maf,this_chromosome):
    #subset on chromosome
    this_index_maf = this_index_maf[this_index_maf.Chromosome == this_chromosome]
    other_maf = other_maf[other_maf.Chromosome == this_chromosome]
    missing_from_index_maf = other_maf[other_maf.Start_Position.isin(this_index_maf.Start_Position) == False]
    return(missing_from_index_maf)


#generate a pileup and get the ref and alt read count for each variant (row) in the missing rows
# add some padding to each region 
# Deletions are handled in a somewhat loosey-goosey way. Basically, the count is the maximum number of reads supporting a deletion at any site
# within the region start-end. This is because there were some examples of a few less reads supporting the deletion near the end (likely the gap in the alignment was shifted)
# DNPs are handled in a super lazy way that should be fine most of the time. The allele is only checked for matching at the first base in the DNP
def get_read_support(chrom,start,end,ref_base,alt_base,bamfile,mut_class,min_base_qual=10, min_mapping_qual=1):
    original_start = start
    actual_start = start - 1 # because of zero vs 1-based indexing
    actual_end = end
    del_positions=[]
    if mut_class == "DEL":
        del_length = len(ref_base)
        del_positions = list(range(actual_start,end))
    if mut_class == "DNP":
        alt_base = alt_base[0]
        ref_base = ref_base[0]
    start = int(start) - padding
    end = int(end) + padding 
    chrom = str(chrom)
    ref_count = 0
    alt_count = 0
    total_depth = 0
    #handle the case in which the bam file is a cram file. This information must be provided by the calling function
    #realname = os.readlink(bamfile)
    realname = str(Path(bamfile).resolve())
    #print(realname)
    if realname.endswith("cram"):
        samfile = pysam.AlignmentFile(realname,'rc')
    else:
        samfile = pysam.AlignmentFile(bamfile,'rb')
    read_bases = {}
    del_support = {}
    max_depth = 0

    for pos in del_positions:
        del_support[pos]=0
    #print(f"{chrom} {start} {end}")
    for pileupcolumn in samfile.pileup(chrom, start, end, ignore_overlaps=False, min_base_quality=min_base_qual,min_mapping_quality=min_mapping_qual,truncate=True):
        if mut_class == "INS":
            #I think this is double-counting overlapping reads but at least it does it consistently
            if pileupcolumn.pos == actual_start:
                seqs = pileupcolumn.get_query_sequences(add_indels=True)
                for base in seqs:
                    if "+" in base:
                        alt_count+=1
                    else:
                        ref_count+=1
                    total_depth+=1
        elif mut_class == "DEL":
            if pileupcolumn.pos in del_positions:
                seqs = pileupcolumn.get_query_sequences(add_indels=True)
                ref_count = 0
                if len(seqs)>max_depth:
                    max_depth = len(seqs)
                for base in seqs:
                    if base == "*":
                        del_support[pileupcolumn.pos]+=1
                    else:
                        ref_count+=1

        elif pileupcolumn.pos == actual_start:
            #print ("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
            for pileupread in pileupcolumn.pileups:
                if mut_class == "SNP" or mut_class == "DNP":
                    if not pileupread.is_del and not pileupread.is_refskip:
                        # query position is None if is_del or is_refskip is set.
                        
                        thisname = pileupread.alignment.query_name
                        this_base = pileupread.alignment.query_sequence[pileupread.query_position]
                        this_qual = pileupread.alignment.query_qualities[pileupread.query_position]
                        if thisname in read_bases:
                            read_bases[thisname].append(this_base)
                        else:
                            read_bases[thisname] = [this_base]
    #now count up bases and deal with disagreement between read 1 and 2 from the same fragment
    if mut_class == "DEL":
        for position in del_support:
            if del_support[position] > alt_count:
                alt_count = del_support[position] 
        ref_count = max_depth - alt_count
        total_depth = max_depth
    if mut_class == "SNP" or mut_class == "DNP":
        for readname in read_bases:
            bases = read_bases[readname]
            #at most there should be two from the same read (fragment/pair etc)
            if len(bases)==1:
                total_depth+=1
                if bases[0] == ref_base:
                    ref_count+=1
                elif bases[0] == alt_base:
                    alt_count+=1
            elif bases[0] == bases[1]:
                total_depth+=1
                if bases[0] == ref_base:
                    ref_count+=1
                elif bases[0] == alt_base:
                    alt_count+=1
            else:
                total_depth+=1
                if bases[0] == ref_base or bases[1] == ref_base:
                    ref_count+=1
                elif bases[0] == alt_base and bases[1] == alt_base:
                    alt_count+=1

        
    #return as strings to prevent pandas from turning them into floats and breaking the MAF in the process
    return({'tumor_depth':str(total_depth),'tumor_ref_count':str(ref_count),'tumor_alt_count':str(alt_count)})

#TODO/WARNING: check if the bam file is a cram file and supply this information to the pileup function. This might break if a cram is provided currently.
# there is code that does this in augment_sv.py but the implementation is more complex here so I haven't gotten to it 
def subset_and_run(full_index_maf,this_chromosome,other_maf_files):
    missing_maf_rows = []
    for other_maf_file in other_maf_files:
        print(f"loading {other_maf_file} for {this_chromosome}")
        other_maf = pd.read_csv(other_maf_file,sep="\t", comment='#',dtype={'Chromosome':'string','t_ref_count':'string','t_alt_count':'string','t_depth':'string'},low_memory=False)
        to_process_maf = get_missing_maf_rows(full_index_maf,other_maf,this_chromosome)
        #check this also against every MAF in missing_maf_rows to avoid processing the same variant more than once
        #print(to_process_maf.shape[0])

        for existing_rows in missing_maf_rows:
            current_dim = to_process_maf.shape[0]
            this_dim = existing_rows.shape[0]
            to_process_maf = get_missing_maf_rows(existing_rows,to_process_maf,this_chromosome)
            #this should remove any overlaping rows 
            new_dim = to_process_maf.shape[0]
            if verbose:
                print(f"subtracting from a maf with {current_dim} using {this_dim} now at {new_dim}")
        missing_maf_rows.append(to_process_maf)

    missing_maf = pd.concat(missing_maf_rows)
    total_to_process = missing_maf.shape[0]
    if verbose:
        print(f"will work on {total_to_process} new MAF rows using {tumour_bam_file} and {normal_bam_file}")
    chroms = missing_maf.Chromosome.tolist()
    start_positions = missing_maf.Start_Position.tolist()
    end_positions = missing_maf.End_Position.tolist()
    ref_alleles = missing_maf.Reference_Allele.tolist()
    alt_alleles = missing_maf.Tumor_Seq_Allele2.tolist()
    mutation_class = missing_maf.Variant_Type.tolist()

    ref_counts = []
    alt_counts = []
    depths = []
    l = len(alt_alleles)
    for i,start_position in enumerate(start_positions):
        read_counts = get_read_support(chrom=chroms[i],
                        start=start_position,
                        end=end_positions[i],
                        ref_base=ref_alleles[i],
                        alt_base=alt_alleles[i],
                        bamfile=tumour_bam_file,
                        mut_class=mutation_class[i])

        ref_counts.append(read_counts['tumor_ref_count'])
        alt_counts.append(read_counts['tumor_alt_count'])
        depths.append(read_counts['tumor_depth'])
        if not i % 500 and verbose:
            if i > 0:
                percent_complete = int(100 * i / l)
                print(f"{i} of {l} rows ({percent_complete}%) complete (chromosome in this thread is {this_chromosome})")
    #add new values to MAF, overwriting what is there
    missing_maf.t_depth = depths
    missing_maf.t_ref_count = ref_counts
    missing_maf.t_alt_count = alt_counts
    return(missing_maf)

def main():
    # load index MAF into a data frame using pandas
    index_maf = pd.read_csv(index_maf_file,sep="\t",   comment='#',dtype={'Chromosome':'string','t_ref_count':'string','t_alt_count':'string','t_depth':'string'},low_memory=False)
    if index_maf_file in add_maf_files:
        index_maf = index_maf.head(1)
    #run in parallel for all chromosomes
    args = []
    for chrom in chromosomes:
        args.append([index_maf,chrom,add_maf_files])

    cool_pool = multiprocessing.Pool(processes = nproc)
    results = cool_pool.starmap(subset_and_run,args)
    cool_pool.close()
    cool_pool.join() #Keeping Chris' variable naming

    #concatenate mafs and write them out to a file
    augmented_maf_merged = pd.concat(results)

    full_augmented_maf_merged = pd.concat([augmented_maf_merged,index_maf])
    full_augmented_maf_merged.Tumor_Sample_Barcode = tumour_sample_id

    #final step: update the Tumor_Sample_Barcode to refer to the index sample
    
    full_augmented_maf_merged.to_csv(augmented_maf,sep="\t",index=False)

main()
