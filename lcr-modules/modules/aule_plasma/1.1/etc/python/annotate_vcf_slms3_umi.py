import argparse 
import pysam

def get_3plus_variants(vcf_3plus):
    """
    Get the variants that are supported by 3 or more variant callers.
    """
    vcf = pysam.VariantFile(vcf_3plus, "r")
    variant_list = []
    for record in vcf:
        # Extract CHR, POS, REF, ALT
        chr = record.chrom
        pos = record.pos
        ref = record.ref
        alt = record.alts[0]  # ALT is a tuple, may contain multiple alternate alleles
        variant_list.append(f"{chr}_{pos}_{ref}_{alt}")
    
    return variant_list

def count_umi_families(bam_file, chrom, pos, ref, alt):
    """
    Count the number of UMI families supporting the variant at the given location.
    Assumes BAM is coordinate-sorted and contains UMI information.
    """
    read_pos = pos - 1  # Convert to 0-based indexing
        
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    umi_family_total = 0
    umi_family_nontrivial = 0
    umi_family_3plus = 0
    umi_family_reference_total = 0
    umi_family_reference_nontrivial = 0
    umi_family_reference_3plus = 0
    umi_family_supporting_total = 0
    umi_family_supporting_nontrivial = 0
    umi_family_supporting_3plus = 0
    
    umi_family_total_fwdstrand = 0
    umi_family_nontrivial_fwdstrand = 0
    umi_family_3plus_fwdstrand = 0
    umi_family_reference_total_fwdstrand = 0
    umi_family_reference_nontrivial_fwdstrand = 0
    umi_family_reference_3plus_fwdstrand = 0
    umi_family_supporting_total_fwdstrand = 0
    umi_family_supporting_nontrivial_fwdstrand = 0
    umi_family_supporting_3plus_fwdstrand = 0
    
    for read in bam.fetch(chrom, read_pos, pos):  # 0-based index
        if read.is_unmapped:
            continue
        
        # get the position in the read where the variant starts
        query_position = [i for i, j in read.get_aligned_pairs() if j == read_pos][0]
        if query_position is None:
            # Position is not covered by this read
            continue
        
        # Find the bases in the read corresponding to the variant
        var_index_in_aligned_pairs = [i for i, j in read.get_aligned_pairs()].index(query_position)
        var_aligned_pairs = read.get_aligned_pairs()[var_index_in_aligned_pairs:(var_index_in_aligned_pairs + max(len(ref), len(alt)))]
        var_seq_indices = [i for i, j in var_aligned_pairs if i is not None]
        base_in_read = "".join([read.query_sequence[i] for i in var_seq_indices]) 
        
        family_size = read.get_tag("cD")
        
        # Count tatal families and separately those that have more than 1 member
        umi_family_total += 1
        if not read.is_reverse:
            umi_family_total_fwdstrand += 1
        if base_in_read == alt:
            umi_family_supporting_total += 1
            if not read.is_reverse:
                umi_family_supporting_total_fwdstrand += 1
        elif base_in_read == ref:
            umi_family_reference_total += 1
            if not read.is_reverse:
                umi_family_reference_total_fwdstrand += 1
        
        if family_size > 1:
            umi_family_nontrivial += 1
            if not read.is_reverse:
                umi_family_nontrivial_fwdstrand += 1
            if base_in_read == alt:
                umi_family_supporting_nontrivial += 1
                if not read.is_reverse:
                    umi_family_supporting_nontrivial_fwdstrand += 1
            elif base_in_read == ref:
                umi_family_reference_nontrivial += 1
                if not read.is_reverse:
                    umi_family_reference_nontrivial_fwdstrand += 1
        
        if family_size > 2:
            umi_family_3plus += 1
            if not read.is_reverse:
                umi_family_3plus_fwdstrand += 1
            if base_in_read == alt:
                umi_family_supporting_3plus += 1
                if not read.is_reverse:
                    umi_family_supporting_3plus_fwdstrand += 1
            elif base_in_read == ref:
                umi_family_reference_3plus += 1
                if not read.is_reverse:
                    umi_family_reference_3plus_fwdstrand += 1
        
    bam.close()
    
    return umi_family_total, umi_family_nontrivial, umi_family_3plus, umi_family_reference_total, umi_family_reference_nontrivial, umi_family_reference_3plus, umi_family_supporting_total, umi_family_supporting_nontrivial, umi_family_supporting_3plus, umi_family_total_fwdstrand, umi_family_nontrivial_fwdstrand, umi_family_3plus_fwdstrand, umi_family_reference_total_fwdstrand, umi_family_reference_nontrivial_fwdstrand, umi_family_reference_3plus_fwdstrand, umi_family_supporting_total_fwdstrand, umi_family_supporting_nontrivial_fwdstrand, umi_family_supporting_3plus_fwdstrand

def annotate_vcf(input_vcf_path, bam_file, vcf_3plus, output_vcf_path):
    """
    Parse the VCF file to extract variant information.
    """
    variants_3plus = get_3plus_variants(vcf_3plus)
    
    vcf_in = pysam.VariantFile(input_vcf_path, "r")
    
    vcf_in.header.info.add("umi_family_total", 1, "Float", "Number of UMI families of any size")
    vcf_in.header.info.add("umi_family_nontrivial", 1, "Float", "Number of UMI families of size 2 or more")
    vcf_in.header.info.add("umi_family_3plus", 1, "Float", "Number of UMI families of size 3 or more")
    vcf_in.header.info.add("umi_family_reference_total", 1, "Float", "Number of UMI families of any size supporting the reference allele")
    vcf_in.header.info.add("umi_family_reference_nontrivial", 1, "Float", "Number of UMI families of size 2 or more supporting the reference allele")
    vcf_in.header.info.add("umi_family_reference_3plus", 1, "Float", "Number of UMI families of size 3 or more supporting the reference allele")
    vcf_in.header.info.add("umi_family_supporting_total", 1, "Float", "Number of UMI families of any size supporting the alternate allele")
    vcf_in.header.info.add("umi_family_supporting_nontrivial", 1, "Float", "Number of UMI families of size 2 or more supporting the alternate allele")
    vcf_in.header.info.add("umi_family_supporting_3plus", 1, "Float", "Number of UMI families of size 3 or more supporting the alternate allele")
    
    vcf_in.header.info.add("umi_family_total_fwdstrand", 1, "Float", "Number of UMI families of any size on forward strand reads")
    vcf_in.header.info.add("umi_family_nontrivial_fwdstrand", 1, "Float", "Number of UMI families of size 2 or more on forward strand reads")
    vcf_in.header.info.add("umi_family_3plus_fwdstrand", 1, "Float", "Number of UMI families of size 3 or more on forward strand reads")
    vcf_in.header.info.add("umi_family_reference_total_fwdstrand", 1, "Float", "Number of UMI families of any size supporting the reference allele on forward strand reads")
    vcf_in.header.info.add("umi_family_reference_nontrivial_fwdstrand", 1, "Float", "Number of UMI families of size 2 or more supporting the reference allele on forward strand reads")
    vcf_in.header.info.add("umi_family_reference_3plus_fwdstrand", 1, "Float", "Number of UMI families of size 3 or more supporting the reference allele on forward strand reads")
    vcf_in.header.info.add("umi_family_supporting_total_fwdstrand", 1, "Float", "Number of UMI families of any size supporting the alternate allele on forward strand reads")
    vcf_in.header.info.add("umi_family_supporting_nontrivial_fwdstrand", 1, "Float", "Number of UMI families of size 2 or more supporting the alternate allele on forward strand reads")
    vcf_in.header.info.add("umi_family_supporting_3plus_fwdstrand", 1, "Float", "Number of UMI families of size 3 or more supporting the alternate allele on forward strand reads")
    
    vcf_in.header.info.add("variant_in_3plus", 0, "Flag", "Whether or not the variant is supported by 3 or more variant callers")
    
    vcf_out = pysam.VariantFile(output_vcf_path, "w", header=vcf_in.header)
    
    for record in vcf_in:
        (
            umi_family_total, umi_family_nontrivial, umi_family_3plus,
            umi_family_reference_total, umi_family_reference_nontrivial, umi_family_reference_3plus,
            umi_family_supporting_total, umi_family_supporting_nontrivial, umi_family_supporting_3plus,
            umi_family_total_fwdstrand, umi_family_nontrivial_fwdstrand, umi_family_3plus_fwdstrand,
            umi_family_reference_total_fwdstrand, umi_family_reference_nontrivial_fwdstrand, umi_family_reference_3plus_fwdstrand,
            umi_family_supporting_total_fwdstrand, umi_family_supporting_nontrivial_fwdstrand, umi_family_supporting_3plus_fwdstrand
        ) = count_umi_families(
            bam_file = bam_file,
            chrom = record.chrom, # Chromosome
            pos = record.pos,     # Position
            ref = record.ref,     # Reference allele
            alt = record.alts[0]  # (first) Alternate allele
        )
        
        variant_in_3plus = record.chrom + "_" + str(record.pos) + "_" + record.ref + "_" + record.alts[0] in variants_3plus
        
        # Add new fields to the INFO dictionary of the record
        record.info["umi_family_total"] = umi_family_total
        record.info["umi_family_nontrivial"] = umi_family_nontrivial
        record.info["umi_family_3plus"] = umi_family_3plus
        record.info["umi_family_reference_total"] = umi_family_reference_total
        record.info["umi_family_reference_nontrivial"] = umi_family_reference_nontrivial
        record.info["umi_family_reference_3plus"] = umi_family_reference_3plus
        record.info["umi_family_supporting_total"] = umi_family_supporting_total
        record.info["umi_family_supporting_nontrivial"] = umi_family_supporting_nontrivial
        record.info["umi_family_supporting_3plus"] = umi_family_supporting_3plus
        
        record.info["umi_family_total_fwdstrand"] = umi_family_total_fwdstrand
        record.info["umi_family_nontrivial_fwdstrand"] = umi_family_nontrivial_fwdstrand
        record.info["umi_family_3plus_fwdstrand"] = umi_family_3plus_fwdstrand
        record.info["umi_family_reference_total_fwdstrand"] = umi_family_reference_total_fwdstrand
        record.info["umi_family_reference_nontrivial_fwdstrand"] = umi_family_reference_nontrivial_fwdstrand
        record.info["umi_family_reference_3plus_fwdstrand"] = umi_family_reference_3plus_fwdstrand
        record.info["umi_family_supporting_total_fwdstrand"] = umi_family_supporting_total_fwdstrand
        record.info["umi_family_supporting_nontrivial_fwdstrand"] = umi_family_supporting_nontrivial_fwdstrand
        record.info["umi_family_supporting_3plus_fwdstrand"] = umi_family_supporting_3plus_fwdstrand
        
        record.info["variant_in_3plus"] = variant_in_3plus
        
        # Write the modified record to the output VCF
        vcf_out.write(record)
    
    # Close the input and output files
    vcf_in.close()
    vcf_out.close()


def main():
    # Parse arguments 
    parser = argparse.ArgumentParser()
    parser.add_argument('-input_vcf_path', metavar='<input_vcf_path>', dest="input_vcf_path", help="input vcf file")
    parser.add_argument('-bam_file', metavar='<bam_file>', dest="bam_file", help="input bam file")
    parser.add_argument('-vcf_3plus', metavar='<vcf_3plus>', dest="vcf_3plus", help="vcf file with variants called by 3 or more variant callers")
    parser.add_argument('-output_vcf_path', metavar='<output_vcf_path>', dest="output_vcf_path", help="output vcf file")
    
    args = parser.parse_args()
  
    # Iterate through the file saving variants metrics to the output file
    annotate_vcf(args.input_vcf_path, args.bam_file, args.vcf_3plus, args.output_vcf_path)

if __name__ == "__main__":
	main()

