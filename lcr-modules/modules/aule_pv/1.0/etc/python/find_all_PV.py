'''
Script to find potential phased variants in a given region and estimate local error rate
'''

import argparse 
import pysam
import csv
import glob
import re
from Bio.Seq import Seq
import pybedtools
from itertools import combinations

import utilis as ut

def merge_sequences(seq1, seq2, overlap):
    if overlap > 0:
        # If there's an overlap, merge by aligning the overlapping regions
        return seq1 + seq2[overlap:]
    else:
        return seq1 + seq2
  
# Function to iterate through provided bed file and parsing statistics
def get_all_pv(bam, min_mapping_qual, min_base_qual, output_doublets, output_triplets, stat, chromsizes):

    # Open BAM file, read all reads
    print("Processing ", bam)
    alignment = pysam.AlignmentFile(bam, "rb")
    reads = alignment.fetch(until_eof=True)

    ## Doublets
    doub_output = output_doublets
    doub_out = open(doub_output, 'w')
    doub_writer = csv.writer(doub_out, delimiter="\t")
    doub_writer.writerow(["PV_ID", "seqnames", "start_pv1", "start_pv2", "ref_pv1", "ref_pv2", "alt_pv1", "alt_pv2", "dist_pv1", "dist_pv2","bq_r1_pv1", "bq_r1_pv2", "bq_r2_pv1", "bq_r2_pv2","both_pv1", "both_pv2","concordance_pv1", "concordance_pv2","rid", "fl", "mq_r1", "mq_r2", "rev_r1", "rev_r2", "r1_all", "r2_all", "r1_n_filtered", "r2_n_filtered", "r1_fs", "r2_fs", "prime5", "prime3", "tl", "r1_len", "r2_len"])

    ## Triplets
    trip_output = output_triplets
    trip_out = open(trip_output, 'w')
    trip_writer = csv.writer(trip_out, delimiter="\t")
    trip_writer.writerow(["PV_ID", "seqnames", "start_pv1", "start_pv2", "start_pv3", "ref_pv1", "ref_pv2", "ref_pv3", "alt_pv1", "alt_pv2", "alt_pv3", "dist_pv1", "dist_pv2", "dist_pv3","bq_r1_pv1", "bq_r1_pv2", "bq_r1_pv3", "bq_r2_pv1", "bq_r2_pv2", "bq_r2_pv3", "both_pv1", "both_pv2", "both_pv3", "concordance_pv1", "concordance_pv2", "concordance_pv3", "rid", "fl", "mq_r1", "mq_r2", "rev_r1", "rev_r2", "r1_all", "r2_all", "r1_n_filtered", "r2_n_filtered", "r1_fs", "r2_fs", "prime5", "prime3", "tl", "r1_len", "r2_len"])

    # # A file with genomic coordinates of active regions
    # bed_out = open(bed, 'w')
    # bed_writer = csv.writer(bed_out, delimiter="\t")

    # Collect metrics
    i = 0  # Total number of reads
    p = 0  # Number of proper pairs
    pq = 0 # Number of proper pairs when both reads are above given quality threshold

    # Assign initial value to r2
    r1 = None
    r2 = None

    # Iterate through reads
    for read in reads:
        # Print progress every 1,000,000 reads
        i += 1
        if i % 1000000 == 0:
            print("%s reads processed" % (i,))

        # Check if read comes from proper pair and not supplementary alignment and is mapped
        if read.is_proper_pair and not read.is_supplementary:

            # Assign reads
            if read.is_read1:
                r1 = read
            elif read.is_read2:
                r2 = read

            # Check if proper pair found
            if not r2 is None and not r1 is None:
                if r1.query_name == r2.query_name:
                    p += 2
                    if r1.mapping_quality > min_mapping_qual and r2.mapping_quality > min_mapping_qual and r1.has_tag("MD") and r2.has_tag("MD"):
                        pq += 2

                        # Proceed if there are any non-reference bases in the read based on the MD tag
                        if ut.countNonRef(r1) > 1 or ut.countNonRef(r2) > 1:
                          # Characterise non-reference bases:
                          # chrom, position, REF, ALT, REF_original (reverse complement if on negative strand), ALT_original (reverse complement if on negative strand)', distance from 5' end of an aligned read, base quality, type ('M' - mismatch, 'I' - insertion, 'D' - deletion) c', 'T', 'c', 'T', 9, 34, 'M'])
                            r1_nonref = ut.parseNonRef(r1)
                            r2_nonref = ut.parseNonRef(r2)

                            discordant_positions = ut.find_discordant_positions(r1, r2) # This returns a list of genomic position, eg [32932537, 32932603, 32932639]
                            
                            if r1_nonref is not None:
                                r1_all = len(r1_nonref)
                                r1_filtered = dict(filter(lambda elem: 
                                  elem[1][7]/len(elem[1][3]) > min_base_qual and elem[1][8] == "M" and elem[1][1] not in discordant_positions, r1_nonref.items()))        # Include only mismatches and filter on BQ
                                r1_n_filtered = (len(r1_nonref) - len(r1_filtered)) # Check how many non-ref bases are left
                            else:
                                r1_all = 0
                                r1_filtered = None
                                r1_n_filtered = 0

                            if r2_nonref is not None:
                                r2_all = len(r2_nonref)
                                r2_filtered = dict(filter(lambda elem: elem[1][7]/len(elem[1][3]) > min_base_qual and elem[1][8] == "M" and elem[1][1] not in discordant_positions, r2_nonref.items()))
                                r2_n_filtered = (len(r2_nonref) - len(r2_filtered))
                            else:
                                r2_all = 0
                                r2_filtered = None
                                r2_n_filtered = 0

                            # Proceed if the number of kept nonreference bases > 1
                            if r1_filtered is not None and r2_filtered is not None:
                                if len(r1_filtered) > 1 and len(r2_filtered) > 1:

                                    # Append read level info (read id, R1|R2 fragment length)
                                    fl = ut.fragLen(r1,r2)        # Aligned fragment length
                                    tl = abs(r1.template_length)  # Sequenced fragment length
                                    rid = str(p)                  # Read ID
                                    mq_r1 = r1.mapping_quality    # MAPQ R1
                                    mq_r2 = r2.mapping_quality    # MAPQ R2
                                    rev_r1 = int(r1.is_reverse)   # Is R1 reverse
                                    rev_r2 = int(r2.is_reverse)   # Is R2 reverse
                                    
                                    # Find fragment sequence 
                                    if (r1.reference_start <= r2.reference_start):
                                      overlap = r1.reference_end - r2.reference_start
                                      complete_sequence = merge_sequences(r1.query_sequence, r2.query_sequence, overlap)
                                    elif (r1.reference_start > r2.reference_start):
                                      overlap = r2.reference_end - r1.reference_start
                                      complete_sequence = merge_sequences(r2.query_sequence, r1.query_sequence, overlap)

                                    # Get the most 5'end and 3'end 4-mers
                                    first_four = complete_sequence[:4]
                                    last_four = complete_sequence[-4:]
                                    
                                    # NOT CURRENTLY INCLUDED - Check for soft clipped bases 
                                    # r1_scb = ut.locateSoftClips(r1)
                                    # r2_scb = ut.locateSoftClips(r2)
                                    # if r1_scb is not None:
                                    #   print(r1_scb)
                                    #   print(complete_sequence)
                                    #   print(["R1:", r1.reference_name, r1.reference_start, r1.reference_end])
                                    #   print(["R2:", r2.reference_name, r2.reference_start, r2.reference_end])
                                    #   print(["TL:", r1.template_length, "FL:",fl])
                                       
                                    # Assign directionality, reverse complement if needed
                                    if r1.is_reverse:
                                      prime5 = ut.complement(last_four)[::-1]
                                      prime3 = ut.complement(first_four)[::-1]
                                    else:
                                      prime5 = first_four
                                      prime3 = last_four

                                    # Save family size info
                                    if r1.has_tag("cD") and r2.has_tag("cD"):
                                        r1_fs = r1.get_tag("cD")
                                        r2_fs = r2.get_tag("cD")
                                    else:
                                        r1_fs = None
                                        r2_fs = None
                                        
                                    # Save length of individual reads 
                                    r1_len = r1.query_length
                                    r2_len = r2.query_length

                                    # Combine read level info into a list
                                    readsinfo = [rid, fl, mq_r1, mq_r2, rev_r1, rev_r2, r1_all, r2_all, r1_n_filtered, r2_n_filtered, r1_fs, r2_fs, prime5, prime3, tl, r1_len, r2_len]

                                    # Make combinations from variants found in R1 and R2
                                    bothreads = ut.intersection(r1_filtered.keys(), r2_filtered.keys())
                                    r1_pvs = r1_filtered.keys()
                                    r2_pvs = r2_filtered.keys()

                                    doublets_both = list(combinations(bothreads, 2))
                                    doublets_r1 = list(combinations(r1_pvs, 2))
                                    doublets_r2 = list(combinations(r2_pvs, 2))

                                    triplets_both = list(combinations(bothreads, 3))
                                    triplets_r1 = list(combinations(r1_pvs, 3))
                                    triplets_r2 = list(combinations(r2_pvs, 3))

                                    # Save all doublets
                                    ## Present at both reads
                                    for pv in doublets_both:
                                        
                                        # PV ID
                                        elems = [r1_filtered[pv[0]][0], \
                                                 r1_filtered[pv[0]][1], r1_filtered[pv[1]][1], \
                                                 r1_filtered[pv[0]][3], r1_filtered[pv[1]][3]]
                                        pv_id = "_".join(map(str, elems))

                                        # Compose a row:
                                        # PV_ID, chrom, pv1_loc, pv2_loc,
                                        # ref1, ref2, alt1, alt2,
                                        # dist1, dist2,
                                        # bq_r1_1, bq_r1_2, bq_r2_1, bq_r2_2
                                        # pv1_both, pv2_both, 
                                        # pv1_discordance, pv2_discordance

                                        row = [pv_id, r1_filtered[pv[0]][0], r1_filtered[pv[0]][1], r1_filtered[pv[1]][1], \
                                               r1_filtered[pv[0]][2], r1_filtered[pv[1]][2], r1_filtered[pv[0]][3], r1_filtered[pv[1]][3], \
                                               r1_filtered[pv[0]][6], r1_filtered[pv[1]][6], \
                                               r1_filtered[pv[0]][7], r1_filtered[pv[1]][7], r2_filtered[pv[0]][7], r2_filtered[pv[1]][7], \
                                               1, 1, \
                                               r1_filtered[pv[0]][9], r1_filtered[pv[1]][9]] 
                                               
                                        doub_writer.writerow(row+readsinfo)

                                        # Write genomic coordinates
                                        pv1_coord = [row[1], row[2]-1, row[2]] # First variant in a doublet
                                        pv2_coord = [row[1], row[3]-1, row[3]] # Second variant in a doublet
                                        #bed_writer.writerow(pv1_coord)
                                        #bed_writer.writerow(pv2_coord)

                                    # R1
                                    for pv in doublets_r1:
                                        if pv[0] not in bothreads or pv[1] not in bothreads:
                                           
                                            # PV ID
                                            elems = [r1_filtered[pv[0]][0], \
                                                     r1_filtered[pv[0]][1], r1_filtered[pv[1]][1], \
                                                     r1_filtered[pv[0]][3], r1_filtered[pv[1]][3]]
                                            pv_id = "_".join(map(str, elems))
  
                                            row = [pv_id, r1_filtered[pv[0]][0], r1_filtered[pv[0]][1], r1_filtered[pv[1]][1], \
                                                   r1_filtered[pv[0]][2], r1_filtered[pv[1]][2], r1_filtered[pv[0]][3], r1_filtered[pv[1]][3], \
                                                   r1_filtered[pv[0]][6], r1_filtered[pv[1]][6]] 
  
                                            # Add base qualities (0 - placeholders for non-overlaping parts)
                                            if pv[0] in bothreads and pv[1] not in bothreads:
                                                bq = [r1_filtered[pv[0]][7], r1_filtered[pv[1]][7], r2_filtered[pv[0]][7],0,1,0] 
                                            elif pv[0] not in bothreads and pv[1] in bothreads:
                                                bq = [r1_filtered[pv[0]][7], r1_filtered[pv[1]][7], 0,r2_filtered[pv[1]][7],0,1]
                                            elif pv[0] not in bothreads and pv[1] not in bothreads:
                                                bq = [r1_filtered[pv[0]][7], r1_filtered[pv[1]][7],0,0,0,0]
  
                                            # Add discordance rate
                                            discord = [r1_filtered[pv[0]][9], r1_filtered[pv[1]][9]]
                                            
                                            # Write row
                                            doub_writer.writerow(row+bq+discord+readsinfo)
  
                                            # Write genomic coordinates
                                            pv1_coord = [row[1], row[2]-1, row[2]] # First variant in a doublet
                                            pv2_coord = [row[1], row[3]-1, row[3]] # Second variant in a doublet
                                            #bed_writer.writerow(pv1_coord)
                                            #bed_writer.writerow(pv2_coord)
                                    # R2
                                    for pv in doublets_r2:
                                        if pv[0] not in bothreads or pv[1] not in bothreads:
                                            # PV ID
                                            elems = [r2_filtered[pv[0]][0], \
                                                     r2_filtered[pv[0]][1], r2_filtered[pv[1]][1], \
                                                     r2_filtered[pv[0]][3], r2_filtered[pv[1]][3]]
                                            pv_id = "_".join(map(str, elems))
  
                                            row = [pv_id, r2_filtered[pv[0]][0], r2_filtered[pv[0]][1], r2_filtered[pv[1]][1], \
                                                   r2_filtered[pv[0]][2], r2_filtered[pv[1]][2], r2_filtered[pv[0]][3], r2_filtered[pv[1]][3], \
                                                   r2_filtered[pv[0]][6], r2_filtered[pv[1]][6]]
  
                                            # Add base qualities
                                            if pv[0] in bothreads and pv[1] not in bothreads:
                                                bq = [r1_filtered[pv[0]][7], 0, r2_filtered[pv[0]][7], r2_filtered[pv[1]][7],1,0]
                                            elif pv[0] not in bothreads and pv[1] in bothreads:
                                                bq = [0, r1_filtered[pv[1]][7], r2_filtered[pv[0]][7], r2_filtered[pv[1]][7],0,1]
                                            elif pv[0] not in bothreads and pv[1] not in bothreads:
                                                bq = [0, 0, r2_filtered[pv[0]][7], r2_filtered[pv[1]][7],0,0]
  
                                            # Add discordance rate
                                            discord = [r2_filtered[pv[0]][9], r2_filtered[pv[1]][9]]
                                            
                                            doub_writer.writerow(row+bq+discord+readsinfo)
  
                                            # Write genomic coordinates
                                            pv1_coord = [row[1], row[2]-1, row[2]] # First variant in a doublet
                                            pv2_coord = [row[1], row[3]-1, row[3]] # Second variant in a doublet
                                            #bed_writer.writerow(pv1_coord)
                                            #bed_writer.writerow(pv2_coord)

                                    # Save all triplets
                                    ## Present at both reads
                                    for pv in triplets_both:
                                        # PV ID
                                        elems = [r1_filtered[pv[0]][0], \
                                                 r1_filtered[pv[0]][1], r1_filtered[pv[1]][1], r1_filtered[pv[2]][1], \
                                                 r1_filtered[pv[0]][3], r1_filtered[pv[1]][3], r1_filtered[pv[2]][3]]
                                        pv_id = "_".join(map(str, elems))

                                        # Compose a row:
                                        # PV_ID, chrom, pv1_loc, pv2_loc,  pv3_loc,
                                        # ref1, ref2, ref, alt1, alt2, alt3,
                                        # dist1, dist2, dist3,
                                        # bq_r1_1, bq_r1_2, bq_r1_3, bq_r2_1, bq_r2_2, bq_r2_3
                                        # pv1_both, pv2_both, pv3_both
                                        # pv1_discordance, pv2_discordance, pv3_discordance
                  

                                        row = [pv_id, r1_filtered[pv[0]][0], r1_filtered[pv[0]][1], r1_filtered[pv[1]][1], r1_filtered[pv[2]][1], \
                                                r1_filtered[pv[0]][2], r1_filtered[pv[1]][2], r1_filtered[pv[2]][2], r1_filtered[pv[0]][3], r1_filtered[pv[1]][3], r1_filtered[pv[2]][3], \
                                                r1_filtered[pv[0]][6], r1_filtered[pv[1]][6], r1_filtered[pv[2]][6], \
                                                r1_filtered[pv[0]][7], r1_filtered[pv[1]][7], r1_filtered[pv[2]][7], r2_filtered[pv[0]][7], r2_filtered[pv[1]][7], r2_filtered[pv[2]][7], \
                                                1, 1, 1, \
                                                r1_filtered[pv[0]][9], r1_filtered[pv[1]][9], r1_filtered[pv[2]][9]] 
                                        trip_writer.writerow(row+readsinfo)

                                        # Write genomic coordinates
                                        pv1_coord = [row[1], row[2]-1, row[2]] # First variant in a triplet
                                        pv2_coord = [row[1], row[3]-1, row[3]] # Second variant in a triplet
                                        pv3_coord = [row[1], row[4]-1, row[4]] # Third variant in a triplet
                                        #bed_writer.writerow(pv1_coord)
                                        #bed_writer.writerow(pv2_coord)
                                        #bed_writer.writerow(pv3_coord)

                                    # R1
                                    for pv in triplets_r1:
                                        if pv[0] not in bothreads or pv[1] not in bothreads or pv[2] not in bothreads:
                                            # PV ID
                                            elems = [r1_filtered[pv[0]][0], \
                                                     r1_filtered[pv[0]][1], r1_filtered[pv[1]][1], r1_filtered[pv[2]][1], \
                                                     r1_filtered[pv[0]][3], r1_filtered[pv[1]][3], r1_filtered[pv[2]][3]]
                                            pv_id = "_".join(map(str, elems))

                                            row = [pv_id, r1_filtered[pv[0]][0], r1_filtered[pv[0]][1], r1_filtered[pv[1]][1], r1_filtered[pv[2]][1], \
                                                    r1_filtered[pv[0]][2], r1_filtered[pv[1]][2], r1_filtered[pv[2]][2], r1_filtered[pv[0]][3], r1_filtered[pv[1]][3], r1_filtered[pv[2]][3],\
                                                    r1_filtered[pv[0]][6], r1_filtered[pv[1]][6], r1_filtered[pv[2]][6]]

                                            # Add base qualities
                                            if pv[0] not in bothreads and pv[1] not in bothreads and pv[2] not in bothreads:
                                                bq = [r1_filtered[pv[0]][7], r1_filtered[pv[1]][7], r1_filtered[pv[2]][7], 0, 0, 0, 0, 0, 0]
                                            elif pv[0] not in bothreads and pv[1] not in bothreads and pv[2] in bothreads:
                                                bq = [r1_filtered[pv[0]][7], r1_filtered[pv[1]][7], r1_filtered[pv[2]][7], 0, 0, r2_filtered[pv[2]][7], 0, 0, 1]
                                            elif pv[0] not in bothreads and pv[1] in bothreads and pv[2] not in bothreads:
                                                bq = [r1_filtered[pv[0]][7], r1_filtered[pv[1]][7], r1_filtered[pv[2]][7], 0, r2_filtered[pv[1]][7], 0, 0, 1, 0]
                                            elif pv[0] in bothreads and pv[1] not in bothreads and pv[2] not in bothreads:
                                                bq = [r1_filtered[pv[0]][7], r1_filtered[pv[1]][7], r1_filtered[pv[2]][7], r2_filtered[pv[0]][7], 0, 0, 1, 0, 0]
                                            elif pv[0] in bothreads and pv[1] in bothreads and pv[2] not in bothreads:
                                                bq = [r1_filtered[pv[0]][7], r1_filtered[pv[1]][7], r1_filtered[pv[2]][7], r2_filtered[pv[0]][7], r2_filtered[pv[1]][7], 0, 1, 1, 0]
                                            elif pv[0] not in bothreads and pv[1] in bothreads and pv[2] in bothreads:
                                                bq = [r1_filtered[pv[0]][7], r1_filtered[pv[1]][7], r1_filtered[pv[2]][7], 0, r2_filtered[pv[1]][7], r2_filtered[pv[2]][7], 0, 1, 1]
                                            elif pv[0] in bothreads and pv[1] not in bothreads and pv[2] in bothreads:
                                                bq = [r1_filtered[pv[0]][7], r1_filtered[pv[1]][7], r1_filtered[pv[2]][7], r2_filtered[pv[0]][7], 0, r2_filtered[pv[2]][7], 1, 0, 1]

                                            discord = [r1_filtered[pv[0]][9], r1_filtered[pv[1]][9], r1_filtered[pv[2]][9]]
                                            
                                            trip_writer.writerow(row+bq+discord+readsinfo)

                                            # Write genomic coordinates
                                            pv1_coord = [row[1], row[2]-1, row[2]] # First variant in a triplet
                                            pv2_coord = [row[1], row[3]-1, row[3]] # Second variant in a triplet
                                            pv3_coord = [row[1], row[4]-1, row[4]] # Third variant in a triplet
                                            #bed_writer.writerow(pv1_coord)
                                            #bed_writer.writerow(pv2_coord)
                                            #bed_writer.writerow(pv3_coord)

                                    # R2
                                    for pv in triplets_r2:
                                        if pv[0] not in bothreads or pv[1] not in bothreads or pv[2] not in bothreads:
                                            # PV ID
                                            elems = [r2_filtered[pv[0]][0], \
                                                     r2_filtered[pv[0]][1], r2_filtered[pv[1]][1], r2_filtered[pv[2]][1], \
                                                     r2_filtered[pv[0]][3], r2_filtered[pv[1]][3], r2_filtered[pv[2]][3]]
                                            pv_id = "_".join(map(str, elems))

                                            row = [pv_id, r2_filtered[pv[0]][0], r2_filtered[pv[0]][1], r2_filtered[pv[1]][1], r2_filtered[pv[2]][1], \
                                                    r2_filtered[pv[0]][2], r2_filtered[pv[1]][2], r2_filtered[pv[2]][2], r2_filtered[pv[0]][3], r2_filtered[pv[1]][3], r2_filtered[pv[2]][3],\
                                                    r2_filtered[pv[0]][6], r2_filtered[pv[1]][6], r2_filtered[pv[2]][6]]

                                            # Add base qualities
                                            if pv[0] not in bothreads and pv[1] not in bothreads and pv[2] not in bothreads:
                                                bq = [0, 0, 0, r2_filtered[pv[0]][7], r2_filtered[pv[1]][7], r2_filtered[pv[2]][7], 0, 0, 0]
                                            elif pv[0] not in bothreads and pv[1] not in bothreads and pv[2] in bothreads:
                                                bq = [0, 0, r1_filtered[pv[2]][7], r2_filtered[pv[0]][7], r2_filtered[pv[1]][7], r2_filtered[pv[2]][7], 0, 0, 1]
                                            elif pv[0] not in bothreads and pv[1] in bothreads and pv[2] not in bothreads:
                                                bq = [0, r1_filtered[pv[1]][7], 0, r2_filtered[pv[0]][7], r2_filtered[pv[1]][7], r2_filtered[pv[2]][7], 0, 1, 0]
                                            elif pv[0] in bothreads and pv[1] not in bothreads and pv[2] not in bothreads:
                                                bq = [r1_filtered[pv[0]][7], 0, 0, r2_filtered[pv[0]][7], r2_filtered[pv[1]][7], r2_filtered[pv[2]][7], 1, 0, 0]
                                            elif pv[0] in bothreads and pv[1] in bothreads and pv[2] not in bothreads:
                                                bq = [r1_filtered[pv[0]][7], r1_filtered[pv[1]][7], 0, r2_filtered[pv[0]][7], r2_filtered[pv[1]][7], r2_filtered[pv[2]][7], 1, 1, 0]
                                            elif pv[0] not in bothreads and pv[1] in bothreads and pv[2] in bothreads:
                                                bq = [0, r1_filtered[pv[1]][7], r1_filtered[pv[2]][7], r2_filtered[pv[0]][7], r2_filtered[pv[1]][7], r2_filtered[pv[2]][7], 0, 1, 1]
                                            elif pv[0] in bothreads and pv[1] not in bothreads and pv[2] in bothreads:
                                                bq = [r1_filtered[pv[0]][7], 0, r1_filtered[pv[2]][7], r2_filtered[pv[0]][7], r2_filtered[pv[1]][7], r2_filtered[pv[2]][7], 1, 0, 1]

                                            discord = [r2_filtered[pv[0]][9], r2_filtered[pv[1]][9], r2_filtered[pv[2]][9]]
                                            
                                            trip_writer.writerow(row+bq+discord+readsinfo)

                                            # Write genomic coordinates
                                            pv1_coord = [row[1], row[2]-1, row[2]] # First variant in a triplet
                                            pv2_coord = [row[1], row[3]-1, row[3]] # Second variant in a triplet
                                            pv3_coord = [row[1], row[4]-1, row[4]] # Third variant in a triplet
                                            #bed_writer.writerow(pv1_coord)
                                            #bed_writer.writerow(pv2_coord)
                                            #bed_writer.writerow(pv3_coord)

                r1 = None
                r2 = None

    #bed_out.close()
    doub_out.close()
    trip_out.close()

    # Sort and merge active regions
    #allregions = pybedtools.BedTool(bed)
    #merged = allregions.sort().slop(b=20, g=chromsizes).sort().merge().moveto(bed)
    #merged = allregions.sort().merge().slop(b=20, g=chromsizes).moveto(bed)

    # Save metric
    stat_out = open(stat, 'w')
    stat_writer = csv.writer(stat_out, delimiter="\t")
    stat_writer.writerow(["Metric", "Value"])
    stat_writer.writerow(["Proper_pairs", p])
    stat_writer.writerow(["Proper_pairs_qual", pq])
    stat_out.close()

    return(print("Done! Output files:", doub_output, trip_output, stat))

def main():

    # Parse arguments 
    parser = argparse.ArgumentParser()
    parser.add_argument('-bam', metavar='<in_bam>', dest="bam", help="input bam file")
    parser.add_argument('-chromsizes', metavar='<chromsizes>', dest="chromsizes", help="chrom.sizes file to define the chromosome lengths for a given genome")
    parser.add_argument('-output_doublets', metavar='<output_doublets>', dest="output_doublets", help="output file with all phased variant doublets")
    parser.add_argument('-output_triplets', metavar='<output_triplets>', dest="output_triplets", help="output file with all phased variant triplets")
    parser.add_argument('-stat', metavar='<stat>', dest="stat", help="output statistic")
    #parser.add_argument('-bed', metavar='<bed>', dest="bed", help="bed file with genomic coordinates containing phased varaints")
    parser.add_argument('-min_mapping_qual', metavar='<min_mapping_qual>', dest="min_mapping_qual", type = int, help="minimum reads mapping quality")
    parser.add_argument('-min_base_qual', metavar='<min_base_qual>', dest="min_base_qual", type = int, help="minimum base quality")

    args = parser.parse_args()
  
    # Iterate through the file saving variants metrics to the output file
    get_all_pv(args.bam, args.min_mapping_qual, args.min_base_qual, args.output_doublets, args.output_triplets, args.stat, args.chromsizes)

if __name__ == "__main__":
	main()
