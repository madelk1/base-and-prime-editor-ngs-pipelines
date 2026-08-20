# RNA-Seq Off-Target Analysis: ABE8e Base Editors

Comprehensive transcriptome-wide off-target analysis of adenine base editors.

## Overview

This repository contains the complete GATK4-based pipeline for detecting and characterizing off-target base editing events in RNA-seq data. The analysis identifies high-confidence off-target sites using stringent multi-stage filtering and replicate-based thresholds to distinguish genuine editing events from technical noise.

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

## Results Summary

### RNP vs mRNA (PCSK9 target, HUES64)

| Stage | RNP | mRNA | Ratio |
|-------|-----|------|-------|
| Individual samples | 16,718 | 26,109 | 1.6× |
| Pooled (2/3 reps) | 100 | 375 | 3.8× |
| Combined 3/6 threshold | 34 | 34 | 1.0× |

**Key insight:** Complete non-overlap (0%) indicates delivery method fundamentally determines off-target profile.

### Architectural Comparison (HiFiRNASeq, HEK293T)

| Architecture | High-Confidence Sites | vs ABE8e | Mean Editing |
|---|---|---|---|
| ABE8e | 2,522 | — | 20.6% |
| HiFi-ABE | 263 | 9.6× improvement | 17.9% |
| Inlaid-ABE | 461 | 5.5× improvement | 18.6% |

## Usage

### Run Full Pipeline on HPC

```bash
# 1. Prepare STAR index
sbatch star_index.slurm

# 2. Align fastqs and preprocess BAMs
sbatch GATK4-fastq-ARbam-hg38.slurm

# 3. Quantify coverage
sbatch bam-readcount-CL.slurm

# 4. Variant calling
sbatch haplotypecaller-gatk4.slurm

# 5. Merge VCF + readcounts
sbatch step2-gatk4.slurm

# 6. Filter variants
sbatch step3-gatk4.slurm

# 7. Multi-threshold analysis
Rscript combined_analysis_2of6_3of6_4of6.R

# 8. Overlap comparison
Rscript compare_conditions_overlap.R

# 9. Annotation
Rscript annotate_offtargets.R
```

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

## Authors

Madeleine K. Lapinaite, Alapati Lab, Arizona State University

## Citation

If you use this pipeline, please cite:
- GATK Best Practices for RNA-seq: https://github.com/broadinstitute/gatk
- Aryee/Joung lab pipeline: https://github.com/aryeelab/RNAseq_BE_editing

## License

MIT

## Contact

For questions or issues, please open an issue on this repository.

---

**Last updated:** August 2026
**Pipeline version:** GATK4 v4.6.2.0
