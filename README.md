# base-and-prime-editor-ngs-pipelines

Code for analyzing next-generation sequencing data from CRISPR base editor and prime editor experiments.

## Overview

This repository contains analysis pipelines for processing and analyzing NGS data from:
Base editor on-target and off-target characterization via Amplicon NGS data.
Prime editor on-target and byproduct characterization via Amplicon NGS data.
Base editor transcriptome-wide off-target effect characterization via total rRNA-depleted RNA-Seq data.

## Repository Structure
main/
  copy_and_rename_fastq.sh
  README.md
  
  base-editor-amplicon-ngs/
    - pipeline_script.sh
    - submit_pipeline_array.sh
    - sample_sheet.csv
    - processed_crispresso_batchfile.txt
    - Batch_Analysis_BE.py
  
  prime-editor-amplicon-ngs/
    copy_and_rename_fastq.sh
    submit_pe_crispresso_array.sh
    Batch_Analysis_PE.py
      
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

## Dependencies
- CRISPResso2
- Python 3.8+
- Pandas
- NumPy
- SLURM job scheduler

## Author
Madeleine King
Arizona State University, 2026  
