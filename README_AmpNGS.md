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
