'''
Useful functions
'''

import pysam
import csv
import glob
import re
from Bio.Seq import Seq

# Compute fragment length
def fragLen(r1, r2):
    r1_start = r1.query_alignment_start
    r1_end = r1.query_alignment_end
    r2_start = r2.query_alignment_end
    r2_end = r2.query_alignment_end

    coord = [r1_start, r1_end, r2_start, r2_end]

    return(max(coord)-min(coord))

# Find complement bases
def complement(base):
    try:
        pairs = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
        compl = "".join(pairs[n] for n in base)
    except KeyError:
        newbase = []
        for n in base:
            if not n.lower() in ['a', 'c', 'g', 't']:
                n = "N"
                newbase.append(n)
            else:
                newbase.append(n)
        compl = "".join(pairs[n] for n in newbase)
    return compl

# Get an intersection of two lists
def intersection(list1, list2):
    intersect = [value for value in list1 if value in list2]
    return intersect

# Function to localise soft clipped bases, returns None if no soft clipped bases found
def locateSoftClips(read):
    if "S" in read.cigarstring:
        cigar = read.cigartuples

        # Check if 5'end soft clipped
        if cigar[0][0] == 4:
            end5 = cigar[0][1]
        else:
            end5 = 0

        # Check if 3'end soft clipped
        if cigar[-1][0] == 4:
            end3 = cigar[-1][1]
        else:
            end3 = 0

        # Return a dictionary of soft clipped bases
        return({"end5": end5, "end3": end3})

    else:
        return None

# Function to count the number of insertions from CIGAR string
def countInsertions(read):
    cigar = read.cigartuples
    ins = 0
    for op in cigar:
        if op[0] == 1:
            ins += 1

    return ins

# Function to to count non-reference position from MD tag
def countNonRef(read):
    # Count mismatches and deletions
    if read.has_tag("MD"):
        MD = read.get_tag("MD")
        mdSub = re.sub(r'([\\^]*[ACGT]+)[0]*', ' \\1 ', MD)
        mdSplit = re.split('[ ]+', mdSub)
        mm = list(filter(lambda elem: re.match(r'^[ACGT]', elem), mdSplit))
        dels = list(filter(lambda elem: re.match(r'^[\\^]', elem), mdSplit))

        # Count insertions
        ins = countInsertions(read)

        # Sum nonref
        nonref = len(mm) + len(dels) + ins
        return nonref
    else:
        return 0

# Function to identify and summarise nonreference bases in a read
def parseNonRef(read):
    
    # Count the number of reference bases (MD tag + CIGAR)
    nonref = countNonRef(read)

    # Process only reads with non-reference bases
    if nonref > 0:
        chrom = read.reference_name
        aligned = read.get_aligned_pairs(with_seq = True)
        seq = read.query_sequence
        qual = read.query_qualities
        rev = int(read.is_reverse)

        # Soft clipped bases to exclude
        softclipped = locateSoftClips(read)
        if not softclipped is None:
            omit5 = softclipped.get("end5")
            omit3 = softclipped.get("end3")
        else:
            omit5 = 0
            omit3 = 0

        # Exclude soft clipped bases (analyse only aligned portion of a read)
        aligned = aligned[omit5:-omit3 or None]

        # Measure distance from 5'end, cound backwards for reverse strand
        if read.is_reverse:
            dist = -len(aligned)
        else:
            dist = 0

        # Remember the last genomic coordinate
        last_pos = aligned[0][0]
        last_id = ""
        # Iterate through the alignment
        # Return nonref bases as list:
        # [chrom, position, REF, ALT, base quality, distance form 5'end]
        # If on reverse strand, reverse complement
        bases = []
        ids = []
        last = []
        for pos in aligned: # Example of a pos (70, 105712216, 'A') - position in the read, coordinate, reference base
            # Insertion
            if pos[1] is None and pos[2] is None:
                alt = seq[pos[0]]
                ref = "."
                if read.is_reverse:
                    alt_orignal = complement(seq[pos[0]])
                    ref_orignal = "."
                else:
                    alt_orignal = seq[pos[0]]
                    ref_orignal = "."

                # Update the last element (if long insertion) or add new
                now = "_".join(map(str, [chrom, last_pos, ref, "Ins"]))
                if now == last:
                    cumqual = bases[-1][7] + qual[pos[0]]
                    idd = "_".join(map(str, [chrom, last_pos, ref, bases[-1][3]+alt]))
                    bases[-1] = [chrom, last_pos, ref, bases[-1][3]+alt, ref_orignal, bases[-1][5]+alt_orignal, abs(dist), cumqual, "I"]
                    ids[-1] = idd
                else:
                    nonref = [chrom, last_pos, ref, alt, ref_orignal, alt_orignal, abs(dist), qual[pos[0]], "I"]
                    idd = "_".join(map(str, [chrom, last_pos, ref, alt]))
                    last = now
                    bases.append(nonref)
                    ids.append(idd)
                dist += 1

            # Deletion
            elif pos[0] is None:
                alt = "."
                ref = pos[2]
                if read.is_reverse:
                    ref_orignal = complement(pos[2])
                else:
                    ref_orignal = pos[2]

                nonref = [chrom, pos[1], ref, alt, ref_orignal, alt, abs(dist), 0, "D"]
                idd = "_".join(map(str, [chrom, pos[1], ref, alt]))
                bases.append(nonref)
                ids.append(idd)
                dist += 1
                last_pos = pos[1]

            # Mismatch
            elif not pos[1] is None and not pos[2] is None and pos[2].islower():
                alt = seq[pos[0]]
                ref = pos[2]
                if read.is_reverse:
                    alt_original = complement(seq[pos[0]])
                    ref_original = complement(pos[2])
                else:
                    alt_original = seq[pos[0]]
                    ref_original = pos[2]
                nonref = [chrom, pos[1], ref, alt, ref_original, alt_original, abs(dist), qual[pos[0]], "M"]
                idd = "_".join(map(str, [chrom, pos[1], ref, alt]))
                bases.append(nonref)
                ids.append(idd)
                dist += 1
                last_pos = pos[1]

            # Move on if canonical
            else:
                dist += 1
                last_pos = pos[1]

        # Pack the list of bases into a dictionary
        all_nonref = dict(zip(ids, bases))

        if len(all_nonref) > 0:
            return(all_nonref)
        else:
            return None

# Reconstruct the aligned sequence from the read, considering CIGAR operations.
def reconstruct_alignment(read):

    reconstructed_seq = ''
    seq_pos = 0

    for op, length in read.cigartuples:
        if op == 0:  # Match or Mismatch
            reconstructed_seq += read.query_sequence[seq_pos:seq_pos+length]
            seq_pos += length
        elif op == 1:  # Insertion
            seq_pos += length
        elif op == 2:  # Deletion
            reconstructed_seq += '-' * length
        elif op == 4:  # Soft clipping
            seq_pos += length
    
    return reconstructed_seq

# Function to find dicordant psoitions between r1 and r2        
def find_discordant_positions(r1, r2):
    
    discordant_positions = []
    
    # Ensure both reads are mapped and paired properly
    if not r1.is_proper_pair or not r2.is_proper_pair:
        return discordant_positions

    # Identify the overlapping region
    overlap_start = max(r1.reference_start, r2.reference_start)
    overlap_end = min(r1.reference_end, r2.reference_end)
    
    # Check for deletions
    r1_aligned = reconstruct_alignment(r1)
    r2_aligned = reconstruct_alignment(r2)
    
    # Check for actual overlap
    if overlap_start >= overlap_end:
        return discordant_positions  # No overlap
    
    # Adjust sequences for overlap
    r1_seq = r1_aligned[overlap_start - r1.reference_start: overlap_end - r1.reference_start]
    r2_seq = r2_aligned[overlap_start - r2.reference_start: overlap_end - r2.reference_start]
    
    # if len(r1_seq) != len(r2_seq) or r1_seq[0] != r2_seq[0]:
    #     print(r1.query_sequence)
    #     print(r2.query_sequence)
    #     print(r1_aligned)
    #     print(r2_aligned)
    #     print(r1_seq)
    #     print(r2_seq)
    #     print(r1.cigarstring)
    #     print(r2.cigarstring)
    #     print(["R1 reference:", r1.reference_start, r1.reference_end])
    #     print(["R2 reference:", r2.reference_start, r2.reference_end])
    #     print(["R1:", overlap_start - r1.reference_start, overlap_end - r1.reference_start])
    #     print(["R2:", overlap_start - r2.reference_start, overlap_end - r2.reference_start])
    #     print(["Overlap start:", overlap_start])
    #     print(["Overlap end:", overlap_end])

    # Check each position for discordance
    for i in range(len(r1_seq)):
        if r1_seq[i] != r2_seq[i]:
            discordant_position = overlap_start + i
            discordant_positions.append(discordant_position)

    return discordant_positions
