# Amplicon DNA NGS of Base Editors and Prime Editors

Code for analyzing next-generation sequencing data from CRISPR base editor and prime editor experiments.

## Overview

This repository contains analysis pipelines for processing and analyzing NGS data from:
- Base editor on-target and off-target characterization via Amplicon NGS data.
- Prime editor on-target and byproduct characterization via Amplicon NGS data.

## Repository Structure
main/
- copy_and_rename_fastq.sh
-  README.md
  
base-editor-amplicon-ngs/
- pipeline_script.sh
- submit_pipeline_array.sh
- sample_sheet.csv
- processed_crispresso_batchfile.txt
- Batch_Analysis_BE.py
  
prime-editor-amplicon-ngs/
- copy_and_rename_fastq.sh
- submit_pe_crispresso_array.sh
- Batch_Analysis_PE.py
      
## Key Differences: Base Editors vs. Prime Editors Amplicon NGS Data

### Base Editor Workflow
1. Rename FASTQ files for simplicity (`copy_and_rename_fastq.sh`)
2. Process reads with respective guide and amplicon information using (`sample_sheet.csv`) through alignment pipeline (`pipeline_script.sh`)
3. Generate CRISPResso batch file with sample info, FASTQ paths, amplicon sequences, and guide sequences (`processed_crispresso_batchfile.txt`)
4. Submit CRISPResso array job via SLURM (`submit_crispresso_array.sh`)
5. Analyze output with `Batch_Analysis_BE.py`

### Prime Editor Workflow
1. Rename FASTQ files for simplicity (`copy_and_rename_fastq.sh`)
2. **Skip pipeline** — use raw data directly to preserve scaffold insertions and indel patterns
3. Generate CRISPResso batch file with sample info, FASTQ paths, guide sequences, PBS, RTT, nicking guide, and scaffold sequences (`rawdata_crispresso_batchfile.txt`)
4. Submit CRISPResso array job via SLURM (`submit_pe_crispresso_array.sh`)
5. Analyze output with `Batch_Analysis_PE.py`

## Input Files

### sample_sheet.csv
Contains metadata for each sample:
- Sample name
- FASTQ file paths (read1, read2)
- Amplicon

### processed_crispresso_batchfile.txt
Contains metadata for each sample:
- Sample name
- FASTQ file paths (R1, R2)
- Amplicon
- Guide RNA

### rawdata_crispresso_batchfile.txt
Contains metadata for each sample:
- Sample name
- FASTQ file paths (read1, read2)
- Amplicon
- Guide RNA (protospacer)
- pegRNA scaffold 
- extension_seq	(pegRNA extension	RTT + PBS combined, 5'→3' as written in the pegRNA (RNA orientation))
- rtt_window_size	(RTT quantification window	Length of the RTT only (not PBS) — e.g. 25 for a 25nt RTT)
- pe_override_seq	(Edited reference sequence	The amplicon sequence after the intended edit is installed — must be ALL UPPERCASE to avoid CRISPResso case-sensitivity errors)
- nicking_guide_seq	(PE3 nicking guide	20nt nicking guide spacer — leave blank for PE2 mode (no nicking guide))

## Output
CRISPResso generates:
- Alignment statistics
- Quantification of on-target edits
- Off-target detection

Python analysis scripts produce:
- Summary tables of editing efficiency
- Off-target burden quantification
- Product distribution analysis (PE only)

## Important Notes
- **Base editors**: Aligned reads are used because preprocessing removes noise without affecting on-target characterization
- **Prime editors**: Raw data is used because alignment can collapse scaffold insertions and indels, losing critical information about termination fidelity
- **FASTQ renaming**: Simplifies batch file generation and reduces filename-related errors
  
## Pipeline Script Breakdown: Amplicon Sequencing Deduplication Workflow

This script processes paired-end amplicon sequencing reads through 7 stages to generate deduplicated FASTQs ready for variant analysis.

---

### Stage 1: UMI Extraction (`umi_tools extract`)

```bash
umi_tools extract \
  -I "$R1" \
  --bc-pattern=NNNNNNNNN \
  --read2-in="$R2" \
  --stdout="$P1" \
  --read2-out="$P2"
```
**What it does:** Extracts the 9-bp **Unique Molecular Identifier (UMI)** from the beginning of Read 1

**Why:** Each DNA fragment gets tagged with a unique barcode during library prep to identify and remove PCR duplicates later

**Parameters:**
- `--bc-pattern=NNNNNNNNN` - 9-bp random barcode position (beginning of Read 1)
- `--read2-in` / `--read2-out` - Process paired reads together
- `--stdout` - Write processed R1 to file

**Output:** UMI moved into read header (appended as `_NNNNNNNNN`) for deduplication tracking

---

### Stage 2: Adapter Trimming (`fastp`)

```bash
fastp \
  -i "$P1" -I "$P2" \
  -o "$TRIMMED_R1" -O "$TRIMMED_R2" \
  --adapter_sequence     AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  --adapter_sequence_r2  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --thread "$THREADS" \
  --html "${SampleID}_fastp.html" \
  --json "${SampleID}_fastp.json"
```

**What it does:** Removes **Illumina sequencing adapters** from the 3' ends of reads

**Why:** Adapters are non-biological sequences added during library prep; they interfere with alignment to the amplicon reference

**Adapters used:**
- **Read 1 adapter:** `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` (TruSeq R1)
- **Read 2 adapter:** `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` (TruSeq R2)

**Parameters:**
- `-i` / `-I` - Input paired-end FASTQs
- `-o` / `-O` - Output trimmed FASTQs
- `--thread` - Parallelization

**Output:** 
- Clean trimmed FASTQs (`${SampleID}_trimmed.1.fastq.gz`, `${SampleID}_trimmed.2.fastq.gz`)
- HTML/JSON QC reports for quality assessment

---

### Stage 3: Amplicon Alignment (`bowtie2`)

```bash
# Build per-sample index (one-time, cached)
bowtie2-build -f "$FASTAS_DIR/${SampleID}.fasta" "$IDX_PREFIX"

# Align reads to amplicon reference
bowtie2 -x "$IDX_PREFIX" \
  -1 "$TRIMMED_R1" -2 "$TRIMMED_R2" \
  -S "$SAM" \
  -p "$THREADS"
```

**What it does:** Aligns cleaned reads to the **expected amplicon sequence**

**Why:** Bowtie2 is optimized for short reads; ensures off-target reads are filtered out

**Reference:** Each sample uses its own amplicon FASTA file (generated from `batch15sample_sheet.csv`)

**Parameters:**
- `-x` - Index prefix (bowtie2 index files)
- `-1` / `-2` - Paired-end reads
- `-S` - Output SAM file
- `-p` - Thread count for alignment

**Per-sample indexing:** Unique index per sample (`${SampleID}_idx_${SLURM_ARRAY_TASK_ID}`) prevents naming conflicts in parallel jobs

**Output:** SAM file (text-based alignment format, ~100-300 MB per sample)

---

### Stage 4: BAM Conversion & Sorting (`samtools`)

```bash
# Convert SAM → BAM (binary, ~6× smaller)
samtools view -@ "$THREADS" -bS "$SAM" -o "$BAM"

# Sort by genomic position (required for deduplication)
samtools sort -@ "$THREADS" "$BAM" -o "$SORTED_BAM"

# Create index for fast random access
samtools index "$SORTED_BAM"
```

**What it does:** Converts text alignment to binary format, sorts, and indexes

**Why:**
- **SAM → BAM:** Reduces file size 4-6×, same information content
- **Sorting:** Deduplication requires reads to be sorted by position
- **Indexing:** Creates `.bai` index for fast queries

**Parameters:**
- `-@` - Thread count for compression/sorting
- `-b` - Output BAM format
- `-S` - Input is SAM

**Output:** Indexed BAM file ready for deduplication

---

### Stage 5: PCR Duplicate Removal (`umi_tools dedup`)

```bash
umi_tools dedup \
  --paired \
  --method=unique \
  --stdin="$SORTED_BAM" \
  --stdout="$DEDUP_BAM"
```

**What it does:** **Removes PCR duplicates** using UMI and alignment position

**How it works:**
1. Groups reads by same position + same UMI tag
2. Keeps only 1 copy per unique (position, UMI) pair
3. Discards duplicates (PCR artifacts)

**Why:** PCR amplification during library prep artificially inflates variant frequencies; UMI-based deduplication recovers true biological signal

**Parameters:**
- `--paired` - Processes paired-end reads
- `--method=unique` - Keep only one read per UMI group


**Result:** Only unique biological molecules retained; PCR bias eliminated

---

### Stage 6: Convert Back to FASTQ (`samtools fastq`)

```bash
# Re-sort by read name (required for paired-end output)
samtools sort -n -@ "$THREADS" \
  -o "${DEDUP_BAM%.bam}.name.bam" "$DEDUP_BAM"

# Extract reads back to paired-end FASTQ
samtools fastq -@ "$THREADS" "${DEDUP_BAM%.bam}.name.bam" \
  -1 "$FINAL_R1" -2 "$FINAL_R2" \
  -0 /dev/null -s /dev/null -n
```

**What it does:** Converts deduplicated BAM back to paired-end FASTQ format

**Why:** Downstream analysis (CRISPResso, variant calling) requires FASTQ input

**Parameters:**
- `-n` - Sort by read name (keeps pairs together)
- `-1` / `-2` - Output paired-end FASTQs
- `-0 /dev/null` - Discard singleton reads

**Output:** Clean deduplicated FASTQs in `dedup_fastq/` directory, ready for analysis

---

### Stage 7: Cleanup

```bash
# Remove intermediate files
rm -f "$SAM" "$BAM" "$P1" "$P2" "$TRIMMED_R1" "$TRIMMED_R2"
rm -f "${IDX_PREFIX}".*.bt2
```

**What it does:** Deletes temporary intermediate files

**Why:** Saves disk space; intermediate files are no longer needed once deduplication completes

---

## Processing Pipeline Summary

| Stage | Tool | Input | Output | Key Function |
|-------|------|-------|--------|--------------|
| 1 | `umi_tools extract` | Raw FASTQ | FASTQ (UMI in header) | Extract 9-bp molecular barcodes |
| 2 | `fastp` | FASTQ | Trimmed FASTQ + QC | Remove Illumina adapters |
| 3 | `bowtie2` | Trimmed FASTQ + Amplicon FASTA | SAM | Align reads to amplicon reference |
| 4 | `samtools` | SAM | Indexed BAM | Convert to binary, sort, index |
| 5 | `umi_tools dedup` | Sorted BAM | Deduplicated BAM | Remove PCR duplicates by UMI |
| 6 | `samtools fastq` | Deduplicated BAM | Final FASTQ | Convert back to FASTQ format |
| 7 | `rm` | Intermediates | — | Clean up temporary files |

## `CRISPResso2 Base Editing`

SLURM array job script for running CRISPResso2 in base editing mode across all samples in a batch.

---

### Overview

Submits one CRISPResso2 job per sample as a SLURM array task. Each task reads its parameters from a row in the batch TSV file and runs CRISPResso2 with base editor output enabled. Unlike prime editing, **aligned and deduplicated reads** are used — the full UMI/fastp/bowtie2/dedup pipeline is run first.

---

### Requirements

- CRISPResso2 v2.3.3 installed in conda environment `ampseq_pipeline`
- Processed (deduplicated) FASTQ files present in the submission directory
- Batch TSV file present in the submission directory
- `logs/` directory (created automatically)

---

### Batch TSV File Format

The script reads a tab-separated file (`processed_crispresso_batchfile.txt`) with the following columns. **One row = one sample.**

| # | Column | Description | Notes |
|---|---|---|---|
| 1 | `name` | Sample label | Spaces OK |
| 2 | `fastq_r1` | R1 FASTQ filename | Deduplicated/processed file |
| 3 | `fastq_r2` | R2 FASTQ filename | Deduplicated/processed file |
| 4 | `amplicon` | Amplicon sequence | Unedited reference |
| 5 | `guide` | Spacer sequence | 20 nt protospacer, no PAM |

---

### CRISPResso2 Parameters

| Parameter | Value | Reason |
|---|---|---|
| `--base_editor_output` | flag | Enables base editing quantification mode |
| `--quantification_window_size` | `20` | Window centered on guide for A-to-G quantification |
| `--quantification_window_center` | `-10` | Centers window at expected edit position (A4–A8 of protospacer) |
| `--exclude_bp_from_left` | `0` | No edge exclusion |
| `--exclude_bp_from_right` | `0` | No edge exclusion |
| `--min_bp_quality_or_N` | `0` | Retains all reads |

---

### Before Submitting

Update these values in the script before each batch:

```bash
#SBATCH --array=1-90        # change to total number of data rows in TSV
TSV="processed_crispresso_batchfile.tsv"   # update to your batch filename
--output_folder "CRISPRessoBatchXX_BEsamples"   # update to batch-specific folder name
```

---

### Submitting

```bash
# From the directory containing processed FASTQs and TSV
sbatch submit_crispresso_array.sh

# Monitor progress
squeue -u madelk1

# Check for errors after completion
grep -r "ERROR\|Traceback" logs/
```

---

### Common Errors

| Error | Cause | Fix |
|---|---|---|
| `⏭️ Skipping` | Row matching failed | Check TSV filename and array range |
| Low read counts | Dedup pipeline discarded too many reads | Check fastp QC report and duplication rate |
| All samples show 0% editing | Wrong amplicon or guide in TSV | Verify sequences against primer design |

---

## `CRISPResso2 Prime Editing `

SLURM array job script for running CRISPResso2 in prime editing mode across all samples in a batch.

---

### Overview

Submits one CRISPResso2 job per sample as a SLURM array task. Each task reads its parameters from a row in the batch TSV file and runs CRISPResso2 with prime editing mode enabled. Raw FASTQ files are used directly — no preprocessing.

---

### Requirements

- CRISPResso2 v2.3.3 installed in conda environment `ampseq_pipeline`
- Raw FASTQ files present in the submission directory
- Batch TSV file present in the submission directory
- `logs/` directory (created automatically)

---

### Batch TSV File Format

The script reads a tab-separated file with the following 10 columns. **One row = one sample.**

| # | Column | Description | Notes |
|---|---|---|---|
| 1 | `name` | Sample label | Spaces OK — CRISPResso converts to underscores |
| 2 | `fastq_r1` | R1 FASTQ filename | Must match exact filename on disk |
| 3 | `fastq_r2` | R2 FASTQ filename | Must match exact filename on disk |
| 4 | `a` | Amplicon sequence | Unedited reference — mixed case OK |
| 5 | `g` | Spacer sequence | 20 nt protospacer, no PAM |
| 6 | `scaffold` | pegRNA scaffold sequence | Leave **blank** for base editor samples — auto-skipped |
| 7 | `extension_seq` | pegRNA 3' extension | Full RTT + PBS, 5'→3' as in pegRNA |
| 8 | `rtt_window_size` | RTT quantification window | Length of **novel sequence only** (not full RTT+PBS) |
| 9 | `pe_override_seq` | Prime-edited reference | Amplicon with edit installed — **must be ALL UPPERCASE** |
| 10 | `nicking_guide_seq` | PE3 nicking guide spacer | Leave **blank** for PE2 mode |

> ⚠️ **`pe_override_seq` must be fully uppercase.** Mixed case causes CRISPResso to fail with `"extension sequence not found in override reference"` due to case-sensitive substring matching.

> ⚠️ **`rtt_window_size` must equal the novel insertion/edit length only**, not the full RTT. For example, a 40 bp loxP insertion uses `rtt_window_size=40` even if the RTT is 74 nt (the remaining 34 nt are a homology arm matching the reference). Setting this incorrectly causes misclassification of all reads as prime-edited.

---

### CRISPResso2 Parameters

| Parameter | Value | Reason |
|---|---|---|
| `--prime_editing_override_prime_edited_ref_seq` | per-sample (col 9) | Provides the expected edited amplicon for alignment |
| `--prime_editing_pegRNA_extension_quantification_window_size` | per-sample (col 8) | Sets quantification window to novel sequence only |
| `--prime_editing_pegRNA_scaffold_min_match_length` | `1` | Detects even minimal scaffold read-through |
| `--min_bp_quality_or_N` | `0` | Retains all reads regardless of base quality |
| `--exclude_bp_from_left` | `0` | No edge exclusion |
| `--exclude_bp_from_right` | `0` | No edge exclusion |
| `--assign_ambiguous_alignments_to_first_reference` | flag | Sends ambiguous reads to unedited reference |
| `--plot_window_size` | `40` | Visualization window |
| `--prime_editing_nicking_guide_seq` | per-sample (col 10) | Added conditionally — only when col 10 is non-empty |

---

### Before Submitting

Update these three values in the script before each batch:

```bash
#SBATCH --array=1-90        # change to total number of data rows in TSV (not counting header)
TSV="rawdata_crispresso_batchfile.tsv"   # update to your batch filename
--output_folder "CRISPRessoBatchXX_PEsamples"   # update to batch-specific folder name
```

---

### Submitting

```bash
# From the directory containing FASTQs and TSV
sbatch submit_pe_crispresso_array.sh

# Monitor progress
squeue -u madelk1

# Check for errors after completion
grep -r "ERROR\|Traceback" logs/
```

---

### Row Matching Logic

Each SLURM array task ID maps to a TSV row by **row number**, not by filename:

```bash
ROW=$(awk -F'\t' -v line="${SLURM_ARRAY_TASK_ID}" 'NR == line+1' "$TSV")
```

Task ID 1 → TSV row 2 (first data row), task ID 2 → TSV row 3, etc. The TSV row order determines which sample each task processes. Ensure rows are ordered correctly before submitting.

---

### Built-in Validation

The script checks the following before running CRISPResso2 and exits with an error message if any check fails:

- `scaffold` column non-empty (skips base editor samples automatically)
- `extension_seq` non-empty
- `rtt_window_size` non-empty
- `pe_override_seq` non-empty
- R1 and R2 filenames reference the same sample number

---

### Common Errors

| Error | Cause | Fix |
|---|---|---|
| `⏭️ Skipping — no scaffold sequence` | Row matching failed — `$ROW` is empty | Check TSV filename in script matches file on disk; check array range |
| `❌ ERROR: extension_seq is empty` | Column 7 blank in TSV | Fill in extension_seq for that sample |
| `"extension sequence not found in override reference"` | `pe_override_seq` contains lowercase | Uppercase the entire override sequence |
| `"calculated prime-edited amplicon is same as reference"` | `pe_override_seq` identical to amplicon | Edit was never added to the override sequence |
| `KeyError: 'U'` | Scaffold contains RNA bases (`U`) | Replace all `U` with `T` in scaffold column |
| `"quantification window out of bounds"` | `--exclude_bp_from_right` overlaps quantification window | Reduce exclusion value or set to 0 |


## `Batch_Analysis_BE.py`

Script for extracting and aggregating base editing outcome metrics across all CRISPResso2 output folders from a single batch run.

---

### Overview

After CRISPResso2 finishes running in base editing mode, this script reads per-sample allele frequency tables, identifies editing events within the protospacer window, and computes standardized outcome metrics across all samples. Unlike the PE analysis script, this script re-analyzes the allele-level data directly rather than reading summary quantification files — it parses each allele sequence character-by-character to classify edits.

---

### Input

For each sample, the script reads:

```
CRISPRessoBatchOutputFolder/
└── CRISPResso_on_<sample_name>/
    └── Alleles_frequency_table_around_sgRNA_<spacer>.txt
```

The sgRNA sequence is extracted directly from the filename. The script automatically searches for this spacer in the reference sequence on both forward and reverse strands to determine editing orientation.

Each row in the allele frequency table represents a unique observed allele with the following columns used:

| Column | Description |
|---|---|
| `Aligned_Sequence` | The aligned allele sequence |
| `Reference_Sequence` | The reference amplicon sequence |
| `#Reads` | Number of reads for this allele |
| `%Reads` | Percentage of total reads |

> Alleles below `CUTOFF_PCT` (default 0.20%) are excluded from all calculations to reduce sequencing noise.

---

### Configuration

Before running, edit **Cell 2** only:

```python
main_folder = "/path/to/CRISPResso/output/folder"

CUTOFF_PCT = 0.20   # minimum allele frequency to include

TARGET_CONFIG = {
    '_Site_3_':  ('reverse', 5),   # reverse strand, target at A5
    '_Site_6_':  ('forward', 7),   # forward strand, target at A7
    '_OT1_':     ('reverse', None) # reverse strand, no specific target (off-target)
}

EDIT_TYPE = "ABE"
SUMMARY_OUTPUT = "editing_analysis_summary.csv"
POSITION_OUTPUT = "editing_positions_for_prism.csv"
```

**`TARGET_CONFIG` keys** must be substrings present in the sample folder name. Each entry specifies:
- **strand**: `'forward'` (A→G editing on + strand) or `'reverse'` (T→C on + strand = A→G on template strand)
- **target_position**: integer A-number from PAM-distal end (e.g. `5` = A5), or `None` for off-target sites with no defined target

> ⚠️ Samples with no matching `TARGET_CONFIG` entry are skipped with a warning. Add all site identifiers before running.

---

### Strand Detection and Position Numbering

The script searches for the spacer sequence in the reference amplicon on both strands. Adenine positions within the protospacer are numbered **PAM-distally (A-numbering)**:

- **Forward strand**: A1 = most PAM-distal adenine (leftmost in sequence as written)
- **Reverse strand**: A1 = most PAM-distal adenine (rightmost in the + strand sequence, since PAM is on the right of the reverse-strand guide)

The configured strand (from `TARGET_CONFIG`) takes precedence over automatic detection.

---

### Output Metrics

All metrics are expressed as **% of alleles ≥ CUTOFF_PCT** (not % of total reads). The script classifies each allele by scanning the aligned sequence within the sgRNA window position-by-position.

| Metric | Definition |
|---|---|
| `Overall_Editing_%` | % of reads containing ≥1 intended edit (A→G on forward, T→C on reverse) anywhere in protospacer |
| `On_Target_%` | % of reads where the designated target adenine (A*N*) is edited, regardless of bystander edits |
| `Purity_%` | % of reads where **only** the target adenine is edited and no other changes exist in the protospacer window (1 difference total) |
| `Unintended_edits_%` | % of reads containing any non-intended substitution (not A→G / T→C) within the protospacer window |

#### Per-position output
A→G (or T→C) editing percentage at each editable position within the protospacer, reported in PAM-distal A-numbering. Separated by strand (forward / reverse) in the Prism output file.

---

### Output Files

| File | Contents |
|---|---|
| `editing_analysis_summary.csv` | Section 1: Overall metrics table (one row per sample). Section 2: Unintended edits table (position, edit type, frequency) |
| `editing_positions_for_prism.csv` | Per-position A→G editing %, pivoted by sample, separated by strand — ready for direct import into GraphPad Prism |

---

### How to Run

```bash
# On HPC — start an interactive session first
srun -A alapinai_lab --pty /bin/bash -i
conda activate ampseq_pipeline

# Run from the directory containing CRISPResso output folders
python Batch_Analysis_BE.py
```

---

### Common Issues

| Symptom | Cause | Fix |
|---|---|---|
| `⚠ No configuration found — SKIPPING` | Sample name doesn't match any `TARGET_CONFIG` key | Add matching pattern to `TARGET_CONFIG` in Cell 2 |
| `On_Target_%` = N/A | `target_position` set to `None` in config | Set target position (e.g. `5` for A5) in `TARGET_CONFIG` |
| `⚠ No allele frequency file found` | CRISPResso job failed or output folder empty | Check SLURM logs; resubmit failed array tasks |
| `⚠ Could not find sgRNA in reference` | Spacer in filename doesn't match amplicon | Verify guide sequence in TSV matches amplicon |
| All samples show 0% editing | Wrong strand configured | Check `TARGET_CONFIG` strand setting |
| Purity always 0% | Multiple edits always present | Check if quantification window is too wide or guide has multiple editable A's |


## `Batch_Analysis_PE.py`

Script for extracting and aggregating prime editing outcome metrics across all CRISPResso2 output folders from a single batch run.

---

### Overview

After CRISPResso2 finishes running in prime editing mode, this script reads the per-sample quantification files and computes standardized editing outcome metrics across all samples. It does **not** re-analyze sequencing reads — it reads the pre-computed tables output by CRISPResso2.

---

### Input

For each sample, the script reads:

```
CRISPRessoBatchOutputFolder/
└── CRISPResso_on_<sample_name>/
    └── CRISPResso_quantification_of_editing_frequency.txt
```

This file is generated automatically by CRISPResso2 and contains a row for each reference category (Reference, Prime-edited, Scaffold-incorporated) with the following columns used by the script:

| Column | Description |
|---|---|
| `Reads_aligned` | Reads assigned to this reference category |
| `Reads_aligned_all_amplicons` | Total reads aligned across all categories (denominator for all % calculations) |
| `Unmodified` | Reads in this category with no additional changes in the quantification window |
| `Modified` | Reads in this category with additional indels or substitutions |
| `Discarded` | Reads filtered by CRISPResso quality thresholds |

---

### Output Metrics

All metrics are expressed as **% of total aligned reads** (`Reads_aligned_all_amplicons`). The script outputs two files:

- **`editing_analysis_summary.csv`** — one row per sample, all metrics as columns
- **`outcomes_for_prism.csv`** — transposed: one row per metric, one column per sample (ready for GraphPad Prism)

#### Metric Definitions

| Metric | Formula | Definition | Reference |
|---|---|---|---|
| `PE_%` | Prime-edited Unmodified / total × 100 | Precise prime editing — intended edit present, no co-occurring indels or substitutions in window. Equivalent to "indel-free editing efficiency" | Doman et al. 2023 |
| `Prime_edits_all_%` | (PE Unmodified + PE Modified) / total × 100 | All reads aligning to prime-edited reference, including those with co-occurring indels. Equivalent to "Prime edits all" | Antoniou et al. 2025 |
| `Scaffold_incorporated_%` | Scaffold Reads_aligned / total × 100 | Reads containing ≥1 pegRNA scaffold-derived nucleotide beyond the RTT, as a fraction of all aligned reads | Doman et al. 2023 |
| `Scaffold_of_PE_%` | Scaffold / (PE + Scaffold) × 100 | Scaffold contamination rate among editing products only | Antoniou et al. 2025 |
| `Indels_on_reference_%` | Reference Modified / total × 100 | Reads on unedited alleles carrying indels in the quantification window — nick occurred but edit was not installed |  |
| `Indels_on_PE_allele_%` | PE Modified / total × 100 | Reads with intended edit plus a co-occurring indel |  |
| `Total_indels_%` | (Ref Modified + PE Modified) / total × 100 | All reads carrying any indel regardless of editing outcome |  |
| `Reference_unmodified_%` | Ref Unmodified / total × 100 | Unedited alleles with no indels — true unedited fraction |  |

> **Note:** `PE_%` is the primary metric used throughout this work for comparing editor performance. `Scaffold_incorporated_%` follows the Doman et al. definition for cross-study comparability. `Scaffold_of_PE_%` (Antoniou definition) is provided for direct comparison with studies reporting scaffold as a fraction of editing products.

---

### How to Run

```bash
# From the directory containing CRISPResso output folders
python Batch_Analysis_PE.py
```

Edit the `main_folder` variable at the top of the script to point to your CRISPResso output directory before running. Results are written to the same directory.

---

### Common Issues

| Error / symptom | Cause | Fix |
|---|---|---|
| Sample showing 0% for all metrics | CRISPResso job failed or output folder missing | Check SLURM logs; resubmit failed array tasks |
| Controls showing ~100% PE% at large insertion sites | `rtt_window_size` set to full RTT length including homology arm | Set `rtt_window_size` to novel insertion length only |
| All samples showing same values | Script reading wrong output folder path | Check `main_folder` variable |
| Missing `CRISPResso_quantification_of_editing_frequency.txt` | CRISPResso crashed mid-run | Check `CRISPResso_RUNNING_LOG.txt` in the sample folder |
---

## Dependencies
- CRISPResso2
- Python 3.8+
- Pandas
- NumPy
- SLURM job scheduler

## Author
Madeleine King
Arizona State University, 2026  
