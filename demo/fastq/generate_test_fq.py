#!/usr/bin/env python3
import random
import gzip

def random_seq(length, letters="ACGT"):
    """Generate a random sequence of a given length."""
    return ''.join(random.choices(letters, k=length))

def random_quality(length, qual_char='I'):
    """Generate a quality string of a given length with a constant quality character."""
    return qual_char * length

def reverse_complement(sequence):
    """Generate the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))

def load_fasta_sequence(fasta_file):
    """Load the first sequence from a FASTA file."""
    sequence = []
    with gzip.open(fasta_file, "rt") as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                sequence.append(line)
    if not sequence:
        raise ValueError(f"No sequence found in {fasta_file}")
    return "".join(sequence)

def get_random_subsequence(sequence, seq_length=150):
    """Extract a random 150bp subsequence from the given sequence."""
    if len(sequence) < seq_length:
        raise ValueError(f"Sequence is shorter than {seq_length}bp")
    start_pos = random.randint(0, len(sequence) - seq_length)
    return sequence[start_pos:start_pos + seq_length]

def main():
    num_reads = 500
    r1_length = 150  # Length of R1 reads
    r2_length = 150  # Length of R2 reads
    umi_length = 8   # Length of UMI reads (8 nt)
    input_fasta = "Homo_sapiens_ENST00000269305_9_sequence.fa.gz"  # Input FASTA file

    # Load the sequence once
    chromosome_seq = load_fasta_sequence(input_fasta)

    with gzip.open("sample_R1.fastq.gz", "wt") as r1_file, \
         gzip.open("sample_R2.fastq.gz", "wt") as r2_file, \
         gzip.open("sample_UMI.fastq.gz", "wt") as umi_file:
        for i in range(1, num_reads + 1):
            read_id = f"read{i}"

            # Get random 150bp subsequence from the loaded sequence for R1
            seq_r1 = get_random_subsequence(chromosome_seq, r1_length)
            qual_r1 = random_quality(r1_length)
            r1_file.write(f"@{read_id}\n{seq_r1}\n+\n{qual_r1}\n")

            # Write R2 entry as reverse complement of R1
            seq_r2 = reverse_complement(seq_r1)
            qual_r2 = random_quality(r2_length)
            r2_file.write(f"@{read_id}\n{seq_r2}\n+\n{qual_r2}\n")

            # Write UMI entry (still random)
            seq_umi = random_seq(umi_length)
            qual_umi = random_quality(umi_length)
            umi_file.write(f"@{read_id}\n{seq_umi}\n+\n{qual_umi}\n")

if __name__ == "__main__":
    main()
