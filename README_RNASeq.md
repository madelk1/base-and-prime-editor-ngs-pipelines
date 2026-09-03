# RNA-Seq Off-Target Analysis: Base Editors

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

## Setup & Reference Files

### Reference Files Verification

**Before running the pipeline, verify all reference files exist and are complete:**

```bash
#!/bin/bash
# Reference Files Verification Script
# Run before starting alignment — confirms all reference files exist and look correct

REF_DIR="/pub/madelk1/rnaseq/refs/hg38"

echo "=========================================="
echo "HG38 Reference Files Verification"
echo "=========================================="

# 1. Reference FASTA
echo ""
echo "1. Reference FASTA:"
ls -lh ${REF_DIR}/GRCh38.primary_assembly.genome.fa 2>/dev/null || echo "❌ MISSING"
grep -c "^>" ${REF_DIR}/GRCh38.primary_assembly.genome.fa 2>/dev/null && echo "✓ Sequences found" || echo "❌ MISSING"

# 2. FASTA Index (.fai)
echo ""
echo "2. FASTA Index (.fai):"
ls -lh ${REF_DIR}/GRCh38.primary_assembly.genome.fa.fai 2>/dev/null || echo "❌ MISSING"

# 3. Sequence Dictionary (.dict)
echo ""
echo "3. Sequence Dictionary (.dict):"
ls -lh ${REF_DIR}/GRCh38.primary_assembly.genome.dict 2>/dev/null || echo "❌ MISSING"

# 4. GTF Annotation
echo ""
echo "4. GTF Annotation:"
ls -lh ${REF_DIR}/gencode.v44.annotation.gtf 2>/dev/null || echo "❌ MISSING"
echo " Genes: $(grep -c '$\tgene\t' ${REF_DIR}/gencode.v44.annotation.gtf 2>/dev/null || echo '0')"
echo " Transcripts: $(grep -c '$\ttranscript\t' ${REF_DIR}/gencode.v44.annotation.gtf 2>/dev/null || echo '0')"

# 5. STAR Index
echo ""
echo "5. STAR Index:"
ls -d ${REF_DIR}/star_index 2>/dev/null || echo "❌ MISSING"
ls ${REF_DIR}/star_index/Genome ${REF_DIR}/star_index/SA ${REF_DIR}/star_index/SAindex 2>/dev/null >/dev/null && echo "✓ Key files present" || echo "❌ Incomplete"

# 6. dbSNP
echo ""
echo "6. dbSNP:"
ls -lh ${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf 2>/dev/null || echo "❌ MISSING"
ls -lh ${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.idx 2>/dev/null || echo "❌ Index MISSING"

# 7. Chromosome Naming Check (must match across all 3)
echo ""
echo "7. Chromosome Naming Check (must match across all 3):"
echo -n " Reference: "
grep "^>" ${REF_DIR}/GRCh38.primary_assembly.genome.fa | head -1 | cut -d' ' -f1
echo -n " GTF: "
grep -v "^#" ${REF_DIR}/gencode.v44.annotation.gtf | head -1 | cut -f1
echo -n " STAR: "
head -1 ${REF_DIR}/star_index/chrName.txt

echo ""
echo "=========================================="
echo "Verification Complete"
echo "=========================================="
```

---

## SLURM Resource Allocation

### Recommended Settings (HPC)

```bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your-email@asu.edu
#SBATCH -A alapinai_lab           # Allocation/account
#SBATCH -J job-name               # Job name
#SBATCH -p standard               # Partition
#SBATCH -t 04:00:00               # Time limit
#SBATCH --cpus-per-task=16        # CPUs per task
#SBATCH --mem=64G                 # Memory
#SBATCH -o logs/job-%j.out        # Output log
#SBATCH -e logs/job-%j.err        # Error log
```

### Per-Step Resource Requirements

| Step | CPU | Memory | Time | Notes |
|------|-----|--------|------|-------|
| STAR index | 16 | 64 GB | 4 hr | One-time; reusable across projects |
| STAR alignment | 16 | 64 GB | 35 min/sample | Scales with sample count |
| SplitNCigarReads | 4 | 32 GB | 2 hr/sample | Memory-intensive; can't parallelize |
| BQSR (Pass 1) | 8 | 32 GB | 90 min/sample | Can run in parallel |
| BQSR (Pass 2, QC) | 8 | 32 GB | 30 min/sample | Light; quick validation |

---

## Expected Inputs & Outputs

### Before Running

- ✓ Paired-end FASTQs at `/pub/madelk1/rnaseq/fastqs/...`
- ✓ STAR index already built
- ✓ Reference files verified (see verification script above)
- ✓ All paths in script match your system

### Expected Outputs

**Alignment Step:**
- `aligned/${sample}Aligned.sortedByCoord.out.bam` — raw STAR alignment

**GATK4 Preprocessing:**
- `forVariantCalling/${sample}.bam` — final analysis-ready BAM
- `metrics/${sample}.rmdup.log`, `.recal_data.table`, `.post_recal_data.table` — QC metrics

**Variant Calling:**
- `variantCalls/${sample}.vcf.gz` — called variants

---

## Critical Fixes Applied: GATK3 → GATK4

| Issue | GATK3 Original | GATK4 Version | Impact |
|-------|---|---|---|
| STAR MAPQ | (not present) | `--outSAMmapqUnique 60` | Moved to alignment stage; removes post-processing |
| RGSM bug | `-RGSM 1` (hardcoded) | `-RGSM ${sample}` | **🔴 Critical:** Fixes VCF multi-sample headers |
| BQSR Pass 2 input | Attempted `-bqsr-recal-file` on temp BAM | Runs on recalibrated BAM from ApplyBQSR | **🔴 Critical:** Validation now works correctly |
| ApplyBQSR syntax | `PrintReads -BQSR` | `ApplyBQSR --bqsr-recal-file` | Tool renamed; parameter structure changed |
| SplitNCigarReads flags | `-rf ReassignOneMappingQuality -U ALLOW_N_CIGAR_READS` | (removed) | GATK4 handles N-cigars natively |

---

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

---

## Post-Processing R Scripts

After per-site editing frequencies are calculated and filtered, two R scripts are used for downstream analysis and annotation.

---

## `Determining Overlapped Off-Target Sites`

Identifies and visualizes overlap between off-target editing sites detected across different experimental conditions (e.g. delivery methods, editor variants, or treatment groups).

### Input

- Per-sample filtered site tables output from the bam-readcount/filtering pipeline (one file per condition)
- Each file contains columns: `chromosome`, `position`, `strand`, `edit_type`, `editing_frequency`, `coverage`

### What it does

1. Loads filtered site tables for each condition
2. Identifies sites present in each condition at the selected replicate threshold (default 3/6)
3. Computes pairwise and multi-way overlaps between conditions
4. Generates Venn diagram or UpSet plot of site overlap
5. Outputs a merged table of all sites with per-condition editing frequencies for comparison

### Configuration

Edit the top of the script before running:

```r
# Paths to filtered site tables (one per condition)
condition_files <- list(
  "ConditionA" = "path/to/conditionA_sites.txt",
  "ConditionB" = "path/to/conditionB_sites.txt"
)

# Replicate threshold used for filtering
min_samples <- 3
```

### Output

| File | Contents |
|---|---|
| `overlap_venn.pdf` | Venn diagram of site overlap across conditions |
| `overlap_merged_table.csv` | All sites with editing frequency per condition — ready for Prism |
| `shared_sites.csv` | Sites present in all conditions (highest-confidence off-targets) |

### How to Run

```r
# From R or RStudio
source("compare_conditions_overlap.R")

# Or from command line
Rscript compare_conditions_overlap.R
```

---

## `Genomic Annotation of Off-Target Sites`

Annotates filtered off-target editing sites with genomic context including gene name, transcript feature (exon/intron/UTR/intergenic), distance to nearest gene, and known SNP overlap.

### Input

- Filtered site table (e.g. `shared_sites.csv` from `compare_conditions_overlap.R` or directly from filtering pipeline)
- Reference genome annotation (GTF/GFF3)
- Optional: dbSNP VCF for known SNP filtering

### What it does

1. Loads filtered sites and converts to GRanges object
2. Overlaps sites with genome annotation to assign feature class (CDS, UTR3, UTR5, intron, intergenic)
3. Assigns nearest gene name and distance for intergenic sites
4. Optionally flags sites overlapping known SNPs in dbSNP
5. Classifies sites by predicted functional impact (coding, splice site, UTR, intronic, intergenic)
6. Outputs annotated table sorted by editing frequency

### Configuration

Edit the top of the script before running:

```r
# Path to filtered sites table
sites_file <- "shared_sites.csv"

# Reference genome annotation (GTF format)
gtf_file <- "path/to/genome.gtf"

# Optional: dbSNP VCF (set to NULL to skip SNP filtering)
dbsnp_vcf <- "path/to/dbsnp.vcf.gz"   # or NULL

# Genome build
genome_build <- "hg38"
```

### Output

| File | Contents |
|---|---|
| `annotated_offtargets.csv` | Full annotated site table with gene, feature, and SNP overlap |
| `offtarget_feature_summary.csv` | Count of sites per feature class (CDS, UTR, intron, intergenic) |
| `offtarget_feature_pie.pdf` | Pie chart of feature class distribution |

### How to Run

```r
# From R or RStudio
source("annotate_offtargets.R")

# Or from command line
Rscript annotate_offtargets.R
```

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

## References

- [GATK Best Practices RNA-seq](https://github.com/broadinstitute/gatk)
- [Aryee/Joung Lab RNAseq_BE_editing](https://github.com/aryeelab/RNAseq_BE_editing)
- [STAR Manual](https://github.com/alexdobin/STAR)
- GENCODE v44 GTF (GRCh38)
- dbSNP 138 VCF
