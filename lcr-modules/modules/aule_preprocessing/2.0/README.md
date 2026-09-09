# aule_preprocessing 2.0

Module for performing alignment and deduplication of raw sequencing data
compatible with downstream SNVs, SVs and phased variants analysis.

## Modification to a standard GATK Best Practices Workflow

- UMIs are read inline from the R1/R2 read sequence itself (`fgbio
  ExtractUmisFromBam`, read structure `5M2S+T 5M2S+T`) and written to the `RX`
  tag, instead of coming from a separate UMI FASTQ file matched by read order.

## Modification to a standard LCR design

- Hard coded absolute paths to the raw FASTQ files (R1 and R2) for safety
  and easy sample tracking.

-------

## Differences from 1.0, and which version to use

1.0 and 2.0 differ in exactly one respect: how the UMI is obtained from the raw
reads. Every other step (adapter marking, alignment, merge, duplicate marking,
BQSR) is identical between the two versions.

| | 1.0 | 2.0 |
|---|---|---|
| UMI source | Separate UMI FASTQ, matched to R1/R2 by read order | Inline in R1/R2 at fixed offsets (read structure `5M2S+T 5M2S+T`) |
| Tool | `fgbio AnnotateBamWithUmis -t RX -q RQ` | `fgbio ExtractUmisFromBam --molecular-index-tags=ZA ZB --single-tag=RX` |
| Sample sheet | Needs `sample_fastq_R1`, `sample_fastq_R2`, `sample_fastq_UMI` | Needs only `sample_fastq_R1`, `sample_fastq_R2` |
| UMI base-quality tag (`RQ`) | Populated | **Not populated** -- `ExtractUmisFromBam` has no equivalent option (see below) |
| Order-desync risk | Real: R1/R2/UMI FASTQs going out of sync silently miscalls the UMI. 1.0's `--fail-fast true` exists specifically to catch this. | Not applicable -- the UMI comes from the same read as the template, so there is nothing to desync. |

**Use 1.0** if your library prep puts the UMI in a separate index read /
FASTQ file (e.g. a dedicated UMI/index read on the sequencer), and you need
`RQ` (UMI base-quality) downstream.

**Use 2.0** if your library prep puts the UMI inline at the start of R1/R2
(a fixed-length UMI + spacer prefix on the template read itself), which is
increasingly the more common design and is what 2.0 assumes via its
`5M2S+T 5M2S+T` read structure. If your kit uses a different UMI/spacer
length, `--read-structure` in `_aule_preprocessing_extract_umis_from_bam`
will need updating for that assay before use -- it is currently hardcoded to
this specific structure, not exposed via `options.extract_umis_from_bam`.

**Do not run library preps designed for one UMI layout through the other
version's module** -- 1.0 against inline UMIs, or 2.0 against a separate
UMI FASTQ, will both silently produce garbage UMI calls with no error, since
neither tool validates that its assumed UMI layout matches the actual data.

### About the missing `RQ` tag in 2.0

`fgbio ExtractUmisFromBam` extracts only the UMI *sequence*; unlike
`AnnotateBamWithUmis -q RQ`, it has no option to also capture the UMI bases'
own sequencing quality. This is a real capability gap in the tool itself, not
a missed flag -- there is nothing to set that restores this in 2.0's config.
2.0's `merge_bam_alignment` options accordingly no longer include
`--ATTRIBUTES_TO_RETAIN RQ` (retaining a tag that's never written would just
be a silent no-op).

Before switching a pipeline from 1.0 to 2.0, check whether anything
downstream actually reads `RQ`'s value (as opposed to just carrying the tag
along) -- as of this writing, `aule_umi_collapse` (1.0 and 2.0) only lists it
under `--ATTRIBUTES_TO_RETAIN`, with no logic keyed on it, so retaining `RQ`
there is itself already a no-op once fed BAMs from this module's 2.0. If you
have (or add) a UMI-quality-aware consensus or filtering step that does read
`RQ`, 2.0 is not a drop-in replacement for 1.0 for that use case without
adding a way to capture UMI-region base quality before `ExtractUmisFromBam`
clips it off the read.

-------

## Workflow

### 1. Soft links to the input FASTQ files: `_aule_preprocessing_input_fastq`

Creates absolute symlinks between paths specified in the sample sheet for R1
and R2 and the `inputs` folder. This is a standard first rule for an
LCR-like workflow. Unlike 1.0, no UMI FASTQ column is read here.

### 2. Converting R1 and R2 FASTQ files to an unmapped BAM: `_aule_preprocessing_fastq_to_unmapped_bam`

Opportunity to assign proper read group and sample metadata to a sequencing file.
Read group ID should differentiate factors contributing to technical batch effects,
i.e. flow cell lane. These info are not stored in a FASTQ file due to the format
constraints.

See [this GATK blog post](https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently) for details.

TODO: Correct the read group meta (does not include flow cell atm, see `bwa_mem`
module from LCR)

### 3. Extract UMIs from the read sequence: `_aule_preprocessing_extract_umis_from_bam`

Reads the UMI inline from the start of each mate using `fgbio
ExtractUmisFromBam` and writes it to the `RX` tag, consuming the UMI and its
adapter/spacer bases out of the read in the process. Replaces 1.0's
`_aule_preprocessing_annotate_bam_with_umis`.

Options

- `--read-structure=5M2S+T 5M2S+T`: 5 bp UMI + 2 bp spacer, then template, for
  each mate. Hardcoded (see "Differences from 1.0" above), not
  `options.extract_umis_from_bam`.
- `--molecular-index-tags=ZA ZB`: per-mate raw UMI tags (fgbio convention).
- `--single-tag=RX`: consolidated UMI tag, matching 1.0's tag name so
  downstream steps (`mark_duplicates --BARCODE_TAG RX`,
  `merge_bam_alignment --ATTRIBUTES_TO_RETAIN RX`) are unchanged.

No UMI base-quality tag is written -- see "About the missing RQ tag" above.

### 4. Mark Illumina adapters: `_aule_preprocessing_mark_illumina_adapters`

This step adds the `XT` tag to a read record to mark the 5' start position of the
specified adapter sequence and produces a metrics file. By default, the tool uses
Illumina adapter sequences. Adapter sequence can be customised with `FIVE_PRIME_ADAPTER`
or `THREE_PRIME_ADAPTER`, [see](https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently).

TODO: Check if we need to customize adapters or stay with the default Illumina's.

### 5. Prepare BWA inputs: `_aule_preprocessing_bwa_input`

This step converts an uBAM file back to interleaved FASTQ (R1 and R2 in one file).
The advantage is that (from a [GATK's blog post](https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently)):

> by specifying `CLIPPING_ATTRIBUTE=XT` and `CLIPPING_ACTION=2`, `SamToFastq`
changes the quality scores of bases marked by XT to two--a rather low score
in the Phred scale. This effectively removes the adapter portion of
sequences from contributing to downstream read alignment and alignment
scoring metrics.

Options

- `--CLIPPING_ATTRIBUTE XT`: uBAM tag specifying adapter position
- `--CLIPPING_ACTION 2`: How to clip; Phred quality of adapter's bases, [see](https://broadinstitute.github.io/picard/command-line-overview.html#SamToFastq)
- `--INTERLEAVE true`: R1 and R2 in one file - one after another
- `--NON_PF true`: all reads included regardless the quality

### 6. BWA MEM alignment: `_aule_preprocessing_bwa_mem`

Reads alignment to the reference genome with BWA-MEM. Self-explanatory.

Options

- `-M`: to flag shorter split hits as secondary (for Picard compatibility), [see](https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently)
- `-p`: R1 and R2 in one file (interleaved paired reads)

### 7. GATK's Merge uBAM and aligned BAM: `_aule_preprocessing_merge_bam_alignment`

This step merges defined information from the unmapped BAM (step 4) with that of
the aligned BAM (step 6) to conserve read data, e.g. original read information
and base sequencing quality scores.

Specifically, aligned BAM generated lacks read group information and certain
tags--the UQ (Phred likelihood of the segment), MC (CIGAR string for mate) and
MQ (mapping quality of mate) tags.

This step also can readjust secondary alignment, [see](https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently).

Options

- `--CREATE_INDEX true`: Whether to create a `.bai` file
- `--ADD_MATE_CIGAR true`: Adds the mate CIGAR tag (MC)
- `--CLIP_ADAPTERS false`: Whether to clip adapters from the reads
- `--CLIP_OVERLAPPING_READS true`: Soft-clips ends so mates do not extend past each other
- `--INCLUDE_SECONDARY_ALIGNMENTS true`: default, self explanatory
- `--MAX_INSERTIONS_OR_DELETIONS -1`: changed to allow any number of insertions or deletions
- `--PRIMARY_ALIGNMENT_STRATEGY MostDistant`: changed from default BestMapq
- `--ATTRIBUTES_TO_RETAIN RX`: Retains a tag with UMI
- `--ATTRIBUTES_TO_RETAIN XS`: Retain a tag with BWA-MEM's suboptimal alignment scores

No `--ATTRIBUTES_TO_RETAIN RQ` -- see "About the missing RQ tag" above.

TODO: I'm not 100% sure how the `--CLIP_OVERLAPPING_READS true` affects variant
calling (esp for ctDNA), but this is a GATK's setting.

### 8. GATK's Fix NM and UQ tags: `_aule_preprocessing_fix_tags`

This tool takes in a coordinate-sorted SAM or BAM and calculates the NM, MD,
and UQ tags by comparing with the reference. This may be needed when
`MergeBamAlignment` was run with SORT_ORDER other than 'coordinate' and thus
could not fix these tags then.

COMMENT: This *may* be a redundant step, but my understanding is that it's either this
or additional sorting before `MergeBamAlignment`.

### 9. GATK's Mark Duplicates: `_aule_preprocessing_mark_duplicates`

This step flags both PCR and optical duplicate reads with a 1024 (0x400) SAM flag.
It's UMI aware thanks to using a `--BARCODE_TAG` option. After duplicate reads are collected,
the tool differentiates the primary and duplicate reads using an algorithm that
ranks reads by the sums of their base-quality scores (default method). By default,
duplicated reads *are not* removed from the file.

Options

- `--BARCODE_TAG RX`: Specify the tag with UMIs
- `--CREATE_INDEX true`: Whether to index the bam file after run

### 10. GATK's Base Quality Recalibrator: `_aule_preprocessing_base_recalibraton`

In a nutshell, it is a data pre-processing step that detects systematic errors
made by the sequencing machine when it estimates the accuracy of each base call.
Read [this post](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-)
for detailed explanation.

TODO: `--known-sites` options in the config file should be included in the reference
workflow

### 11. GATK's Apply recalibrated base qualities: `_aule_preprocessing_apply_bsqr`

This step applies recalibrated base qualities to a bam file.
