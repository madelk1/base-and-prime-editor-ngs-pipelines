# RNA-Seq Off-Target Analysis: ABE8e Base Editors

Comprehensive transcriptome-wide off-target analysis of adenine base editors.

## Overview

This repository contains the complete GATK4-based pipeline for detecting and characterizing off-target base editing events in RNA-seq data. The analysis identifies high-confidence off-target sites using stringent multi-stage filtering and replicate-based thresholds to distinguish genuine editing events from technical noise. It has been adapted from Aryee/Joung lab pipeline: https://github.com/aryeelab/RNAseq_BE_editing to work with GATK4. I have also included additional scripts for high-confidence off-target sites and comparing between datasets (i.e. different variants or delivery formats)


## Pipeline Architecture

### Step 1: GATK4 Preprocessing
Alignment, duplicate marking, split-N-cigar reads, base recalibration

**Files:**
- `download_refs.slurm` - Download of human genome and dbSNP if using
- `star_index.slurm` - STAR index generation
- `GATK4-fastq-ARbam-hg38.slurm` - Initial alignment with STAR v2.7.10a
- `haplotypecaller-gatk4.slurm` - Variant calling

### Step 2: Merge VCF + Coverage
Combine HaplotypeCaller VCF with bam-readcount nucleotide abundances

**Files:**
- `bam-readcount-CL.slurm` - Coverage quantification
- `step2-gatk4.slurm` - VCF/readcount merging

### Step 3: Variant Filtering
Stringent filtering for high-confidence A→G and T→C variants

**Files:**
- `step3-gatk4.slurm` - Per-sample variant filtering

### Step 4: Multi-Threshold Analysis
Replicate-based filtering at 2/6, 3/6, 4/6 thresholds

**Files:**
- `combined_analysis_2of6_3of6_4of6.R` - Threshold optimization

### Step 5: Overlap & Annotation
Compare off-target profiles across conditions; annotate genomic features

**Files:**
- `compare_conditions_overlap.R` - Pairwise and three-way overlaps
- `annotate_offtargets.R` - GTF-based genomic annotation

## File Structure

abe-rnaseq/
- download_refs.slurm
- star_index.slurm
- GATK4-fastq-ARbam-hg38.slurm
- haplotypecaller-gatk4.slurm
- bam-readcount-CL.slurm
- step2-gatk4.slurm
- step3-gatk4.slurm
- gatk4_3stage_pooled_2of3.R
- combined_analysis_2of6_3of6_4of6.R
- annotate_offtargets.R

## Methods & Technical Details

### GATK4 Preprocessing Pipeline

This pipeline implements [GATK Best Practices for RNA-seq variant calling](https://github.com/broadinstitute/gatk), with specific optimizations for base editor off-target detection.

#### Step 1: STAR Alignment

- **Tool:** STAR v2.7.10a, two-pass mode (`--twopassMode Basic`)
- **Rationale:** Two-pass mode identifies novel splice junctions in first pass, then re-aligns reads using both annotated (GENCODE v44) and discovered junctions for improved accuracy across exon boundaries
- **Key parameters:**
  - `--outFilterMultimapNmax 1` - Retains only uniquely-mapping reads (critical for base editing to avoid ambiguous editing calls on multi-mapped positions)
  - `--outSAMmapqUnique 60` - Assigns MAPQ 60 to unique alignments (replaces deprecated GATK3 MAPQ reassignment)
- **Output:** 
  - ~44-70 million unique read pairs per sample
  - 85-95% alignment rates (typical for human RNA-seq)

#### Step 2: AddOrReplaceReadGroups

Adds GATK-required `@RG` tags for downstream processing:
- `RGID` - Read group identifier
- `RGLB` - Library name
- `RGPL` - Platform (Illumina)
- `RGPU` - Platform unit
- `RGSM` - **Sample name** ⚠️ **Critical fix applied:** Original GitHub script hardcoded `-RGSM 1` for all samples, causing incorrect sample identification in VCF headers. Corrected to `-RGSM ${sample}` for proper multi-sample variant calling.

#### Step 3: MarkDuplicates

```bash
--REMOVE_DUPLICATES true
```

**Rationale for RNA-seq:** PCR duplicates artificially inflate variant allele frequencies, critically affecting base editing frequency measurements. Unlike WGS (which marks but retains duplicates for coverage), RNA-seq removes duplicates because highly-expressed genes naturally produce identical reads that are not PCR artifacts.

**Observed metrics:** 10.5-12.6% duplication rates across samples, indicating high library complexity.

#### Step 4: SplitNCigarReads

- Splits reads spanning splice junctions (indicated by 'N' in CIGAR string) into separate exonic alignments
- Hard-clips sequences overhanging introns

**Why essential for RNA-seq:** Standard variant callers assume contiguous genomic DNA; without splitting, junction-spanning reads would be misinterpreted as large deletions, leading to false off-target calls.

**Output:** ~40 million supplementary alignments per sample (normal for RNA-seq, representing split portions of junction-spanning reads)

#### Step 5-6: Base Quality Score Recalibration (BQSR)

**BaseRecalibrator (Pass 1):**
- Builds recalibration model using dbSNP 138 as known variant positions
- Excludes known SNP sites from error model (assumed true variants)
- Uses non-SNP mismatches to calculate empirical quality score distributions
- Accounts for covariates: read position, base context (dinucleotide), machine cycle

**ApplyBQSR:**
- Applies recalibration model to correct quality scores
- ⚠️ **GATK4 syntax:** `--bqsr-recal-file` (not GATK3's `-BQSR`)
- Outputs analysis-ready BAM with improved base quality accuracy

**BaseRecalibrator (Pass 2):**
- Generates post-recalibration QC report for validation
- Recommended for benchmarking BQSR effectiveness

## Variant Calling & Filtering Strategy

### HaplotypeCaller Configuration

```bash
GATK HaplotypeCaller \
  --standard-min-confidence-threshold-for-calling 20.0 \
  --dont-use-soft-clipped-bases
```

Per-base nucleotide abundances quantified via `bam-readcount v1.0.1` for precise editing frequency calculation.

### Three-Stage Per-Sample Filtering

**Stage 1:** Coverage ≥10× in all treated and control samples

**Stage 2:** Restrict to A→G and T→C edits only
- A→G: primary adenine base editor target
- T→C: reverse complement (same biological event)
- Excludes other SNV types (background mutations, sequencing artifacts)

**Stage 3:** ADAR removal
- Remove sites with >1% editing in all control samples
- Eliminates background adenosine deaminase activity
- Results in 0 ADAR sites detected across conditions ✓

### Multi-Threshold Replicate Analysis

High-confidence sites identified across biological replicates using minimum-sample-support thresholds:

| Threshold | Total Sites | A→G | T→C | Mean Editing | Interpretation |
|-----------|-------------|-----|-----|--------------|-----------------|
| 2/6 | 741 | 365 | 376 | 18.1% | Permissive; includes single-modality artifacts |
| **3/6** | **34** | **17** | **17** | **19.7%** | **Optimal; robust reproducible sites** |
| 4/6 | 1 | — | — | 1.8% | Stringent; background-level editing |

**Threshold Selection Rationale:** The 3/6 threshold was selected as optimal because:
1. ✓ Captures reproducible sites detected across delivery methods or replicates
2. ✓ Achieves mean editing frequency (19.7%) substantially above the single 4/6 site (1.8% ≈ background)
3. ✓ Balances sensitivity against false positive rate
4. ✓ Enriches for genuine base editor activity vs. sequencing noise

---

## Known Issues & Corrections

| Issue | Root Cause | Fix |
|-------|-----------|-----|
| Incorrect sample IDs in VCF headers | RGSM hardcoded to "1" | Changed `-RGSM 1` → `-RGSM ${sample}` |
| Invalid GATK4 syntax | GATK3 to GATK4 migration | Changed `-BQSR` → `--bqsr-recal-file` |
| Redundant MAPQ reassignment | GATK3 parameter removed | Moved to STAR upstream (`--outSAMmapqUnique 60`) |

---

### Key Output Files

- `*all_snvs.stringent.txt` - Per-sample high-confidence variants
- `*_offtarget_sites_3stage_2of3.tsv` - Pooled high-confidence sites
- `combined_RNP_mRNA_3of6_offtargets_annotated.tsv` - Final annotated set

## Requirements

### Software
- STAR v2.7.10a
- GATK4 v4.6.2.0
- samtools
- bam-readcount v1.0.1

### R Packages
```r
tidyverse, readr, dplyr, ggplot2, GenomicRanges, rtracklayer
```
]

## References

- [GATK Best Practices RNA-seq](https://github.com/broadinstitute/gatk)
- [Aryee/Joung Lab RNAseq_BE_editing](https://github.com/aryeelab/RNAseq_BE_editing)
- [STAR Manual](https://github.com/alexdobin/STAR)
- GENCODE v44 GTF (GRCh38)
- dbSNP 138 VCF
