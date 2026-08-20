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


## Methods

### Alignment & Preprocessing
- STAR v2.7.10a (2-pass mode, GENCODE v44, GRCh38)
- GATK4 v4.6.2.0 (MarkDuplicates, SplitNCigarReads, BaseRecalibrator/ApplyBQSR)
- dbSNP v138 for BQSR
- PCR duplicate removal

### Variant Calling
- HaplotypeCaller with default settings (--standard-min-confidence-threshold-for-calling 20.0)
- bam-readcount for per-base nucleotide abundances

### Filtering Strategy

**Per-sample (Step 3):**
- ≥10× coverage in treated and control samples
- Removal of sites with >1% editing in controls (ADAR background)
- Restrict to A→G and T→C (base editor signatures)

**Replicate-based (Multi-threshold):**
- 2/6 threshold: Any 2+ samples (permissive)
- 3/6 threshold: Any 3+ samples (optimal balance)
- 4/6 threshold: 4+ samples (stringent)

### Annotation
- GENCODE v44 GTF for genomic feature assignment
- Classification: protein-coding exon, non-coding exon, intronic, intergenic

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

### Reference Files
- GENCODE v44 GTF (GRCh38)
- STAR index (GRCh38)
- dbSNP v138 VCF

