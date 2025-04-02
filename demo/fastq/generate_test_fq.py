#!/usr/bin/env python3
import random
import gzip

def random_seq(length, letters="ACGT"):
    """Generate a random sequence of a given length."""
    return ''.join(random.choices(letters, k=length))

def random_quality(length, qual_char='I'):
    """Generate a quality string of a given length with a constant quality character."""
    return qual_char * length

def main():
    num_reads = 50
    r1_length = 150  # Length of R1 reads
    r2_length = 150  # Length of R2 reads
    umi_length = 8  # Length of UMI reads (8 nt)

    with gzip.open("sample_R1.fastq.gz", "wt") as r1_file, \
         gzip.open("sample_R2.fastq.gz", "wt") as r2_file, \
         gzip.open("sample_UMI.fastq.gz", "wt") as umi_file:
        for i in range(1, num_reads + 1):
            read_id = f"read{i}"

            # Write R1 entry
            seq_r1 = random_seq(r1_length)
            qual_r1 = random_quality(r1_length)
            r1_file.write(f"@{read_id}\n{seq_r1}\n+\n{qual_r1}\n")

            # Write R2 entry
            seq_r2 = random_seq(r2_length)
            qual_r2 = random_quality(r2_length)
            r2_file.write(f"@{read_id}\n{seq_r2}\n+\n{qual_r2}\n")

            # Write UMI entry
            seq_umi = random_seq(umi_length)
            qual_umi = random_quality(umi_length)
            umi_file.write(f"@{read_id}\n{seq_umi}\n+\n{qual_umi}\n")

if __name__ == "__main__":
    main()
