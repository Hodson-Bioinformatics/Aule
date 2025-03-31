# AULE Preprocessing 

Module for performing alignment and deduplication of DIRECT raw sequencing data 
compatible with downstream SNVs, SVs and phased variants analysis. 

#### Modification to a standard GATK Best Practices Workflow 

- Support for UMIs saved in a separate FASTQ file (integrated to `uBAMs` under 
`RX` tag, alongside with their sequencing qualities in `RQ` tag)

#### Modification to a standard LCR design 

- Hard coded absolute paths to the raw FASTQ files (R1, R2 and UMIs) for safety 
and easy sample tracking 

-------

## Workflow

#### 1. Soft links to the input FASTQ files: `_aule_preproc_input_fastq`

Creates absolute symlinks between paths specified in the sample sheet for R1, 
R2 and UMI and the `inputs` folder. This is a standard first rule for an LCR-like
workflow. 


#### 2. Converting R1 and R2 FASTQ files to an unmapped BAM: `_aule_preproc_fastq_to_unmapped_bam`

Opportunity to assign proper read group and sample metadata to a sequencing file. 
Read group ID should differentiate factors contributing to technical batch effects, 
i.e. flow cell lane. These info are not stored in a FASTQ file due to the format
constraints. 

See [this GATK blog post](https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently) for details. 

TODO: Correct the read group meta (does not include flow cell atm, see `bwa_mem` 
module form LCR) 


#### 3. Annotate each read with UMI: `_aule_annotate_bam_with_umis`

This step annotates each read with a appropriate UMI and its sequencing quality. 

##### Options

- `-t RX`: BAM tag to save a UMI's sequence.   
- `-q RQ `: BAM tag to save a UMI's sequencing quality.  
- `--fail-fast true`: Fail immidiately instead of finishing the task for all reads.  

TODO: Consider adding an additional FASTQ validation step to check if corresponding 
reads in R1, R2 and UMIs file are in the same order. (Although there is no such 
problem in the test set, I remember I noticed one or two samples causing the pipeline
to break because of that. I used `fgbio` to fix that, see `validateFASTQ.smk` from
in `Aule_prototype` repository, but I restrained from this here because this was 
a time expensive step). 



#### 4. Mark Illumina adapters: `_aule_mark_illumina_adapters`  

This step adds the `XT` tag to a read record to mark the 5' start position of the 
specified adapter sequence and produces a metrics file. By default, the tool uses 
Illumina adapter sequences. Adapter sequence can be customised with `FIVE_PRIME_ADAPTER`
or `THREE_PRIME_ADAPTER`, [see](https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently).

TODO: Check if we need to customize adapters or stay with the default Illumina's. 


#### 5. Prepare BWA inputs: `_aule_bwa_input`  

This step converts an uBAM file back to interleaved FASTQ (R1 and R2 in one file). 
The advantage is that (from a [GATK's blog post](https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently)):  

> by specifying `CLIPPING_ATTRIBUTE=XT` and `CLIPPING_ACTION=2`, `SamToFastq` changes 
the quality scores of bases marked by XT to two--a rather low score in the Phred scale. 
This effectively removes the adapter portion of sequences from contributing to 
downstream read alignment and alignment scoring metrics.

##### Options

- `CLIPPING_ATTRIBUTE=XT`: uBAM tag specyfying adapter position   
- `CLIPPING_ACTION=2`: How to clip; Phred quality of adapter's bases, [see](https://broadinstitute.github.io/picard/command-line-overview.html#SamToFastq)
- `INTERLEAVE=true`: R1 and R2 in one file - one after another
-  `NON_PF=true` : all reads included regardless the quality

#### 6. BWA MEM alignment: `_aule_bwa_mem`

Reads alignment to the reference genome with BWA-MEM. Self-explanatory.

##### Options

- `-M`: to flag shorter split hits as secondary (for Picard compatibility), [see](https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently) 
- `-p`: R1 and R2 in one file (interleaved paired reads)


#### 7. GATK's Merge uBAM and aligned BAM: `_aule_merge_bam_alignment`

This step merges defined information from the unmapped BAM (step 4) with that of 
the aligned BAM (step 6) to conserve read data, e.g. original read information 
and base sequencing quality scores.

Specifically, aligned BAM generated lacks read group information and certain 
tags--the UQ (Phred likelihood of the segment), MC (CIGAR string for mate) and 
MQ (mapping quality of mate) tags.  

This step also can readjust secodnary alignment, [see](https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently).  

##### Options

- `--CREATE_INDEX true`: Whether to create a `.bai` file
- `--ADD_MATE_CIGAR true`:  Adds the mate CIGAR tag (MC) 
- `--CLIP_ADAPTERS false`:  Whether to clip adapters from the reads  
- `--CLIP_OVERLAPPING_READS true`: Soft-clips ends so mates do not extend past each other 
- `--INCLUDE_SECONDARY_ALIGNMENTS true`: default, self explanatory
- `--MAX_INSERTIONS_OR_DELETIONS -1`: changed to allow any number of insertions or deletions
- `--PRIMARY_ALIGNMENT_STRATEGY MostDistant`: changed from default BestMapq
- `--ATTRIBUTES_TO_RETAIN RX`: Retains a tag with UMI
- `--ATTRIBUTES_TO_RETAIN RQ`: Retains a tag with UMI qualities
- `--ATTRIBUTES_TO_RETAIN XS`: Retain a tag with BWA-MEM's suboptimal alignment scores

TODO: I'm not 100% sure how the `--CLIP_OVERLAPPING_READS true` affects variant calling (esp for ctDNA),
but this is a GATK's setting. 

#### 8. GATK's Fix NM and UQ tags: `_aule_fix_tags`

This tool takes in a coordinate-sorted SAM or BAM and calculates the NM, MD, 
and UQ tags by comparing with the reference. This may be needed when 
`MergeBamAlignment` was run with SORT_ORDER other than 'coordinate' and thus 
could not fix these tags then. 

JAK COMMENT: This *may* be a redundant step, but my understanding is that it's either this 
or additional sorting before `MergeBamAlignment`.

#### 9. GATK's Mark Duplicates: `_aule_mark_duplicates`

This step flags both PCR and optical duplicate reads with a 1024 (0x400) SAM flag. 
It's UMI aware thanks to using a `--BARCODE_TAG` option. After duplicate reads are collected, 
the tool differentiates the primary and duplicate reads using an algorithm that 
ranks reads by the sums of their base-quality scores (default method). By default,
duplicated reads *are not* removed from the file. 

##### Options

- `--BARCODE_TAG RX`: Specify the tag with UMIs 
- `--CREATE_INDEX true`: Whether to index the bam file after run

#### 10. GATK's Base Quality Recalibrator: `_aule_base_recalibraton`

In a nutshell, it is a data pre-processing step that detects systematic errors 
made by the sequencing machine when it estimates the accuracy of each base call.
Read [this post](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-) 
for detailed explanation. 

TODO: `--known-sites` options in the config file should be included in the reference 
workflow 

#### 11. GATK's Apply recalibrated base qualities: `_aule_apply_bsqr`

This step applies recalibrated base qualities to a bam file. 
