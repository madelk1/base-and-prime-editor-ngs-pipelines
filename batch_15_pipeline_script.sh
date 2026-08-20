#!/bin/bash

set -euo pipefail

# ---------------------------
# Range filter (edit or pass "START-END" as $1, e.g. 211-318)
START=1
END=519
if [[ $# -gt 0 && "$1" =~ ^([0-9]+)-([0-9]+)$ ]]; then
  START="${BASH_REMATCH[1]}"
  END="${BASH_REMATCH[2]}"
fi
echo "Processing SampleIDs in range: ${START}-${END}"
# ---------------------------

# Threads
THREADS="${SLURM_CPUS_PER_TASK:-4}"

# Load environments
module load mamba/latest
source activate ampseq_pipeline

# File paths
CSV="batch15sample_sheet.csv"
FASTAS_DIR="amplicon_fastas"
OUTPUT_FASTQ_DIR="dedup_fastq"

mkdir -p "$FASTAS_DIR" "$OUTPUT_FASTQ_DIR" logs

# Helper to trim CR characters (in case CSV has Windows line endings)
trim_cr() { sed $'s/\r$//' ; }

# ===== Generate FASTA files (only for rows in range) =====
echo "📄 Generating FASTA files using SampleID..."
tail -n +2 "$CSV" | trim_cr | while IFS=',' read -r SampleID Index Guide Amplicon; do
  # Basic checks
  if [[ -z "${SampleID:-}" || -z "${Amplicon:-}" ]]; then
    echo "Skipping due to missing SampleID or Amplicon: $SampleID"
    continue
  fi

  # Enforce numeric range
  if [[ "$SampleID" =~ ^[0-9]+$ ]]; then
    if (( SampleID < START || SampleID > END )); then
      echo "⏭️  Skip $SampleID (outside ${START}-${END})"
      continue
    fi
  else
    # Non-numeric IDs are skipped when a range is active
    echo "⏭️  Skip non-numeric SampleID '$SampleID' due to active numeric range"
    continue
  fi

  fasta_file="$FASTAS_DIR/${SampleID}.fasta"
  # (Re)create FASTA file
  echo ">${SampleID}" > "$fasta_file"

  # Clean and validate amplicon sequence
  clean_amplicon=$(echo "$Amplicon" | tr -d '\r\n[:space:]' | tr 'acgt' 'ACGT' | tr -cd 'ACGT')
  if [[ -z "$clean_amplicon" ]]; then
    echo "ERROR: Empty or invalid amplicon for $SampleID"
    rm -f "$fasta_file"
    continue
  fi
  echo "$clean_amplicon" >> "$fasta_file"
done

# ===== Bowtie2 indexes will be built per-sample with a unique prefix =====
echo "🔧 Using per-sample Bowtie2 indexes..."

# ===== Run processing =====
echo "🚀 Running processing on samples..."
tail -n +2 "$CSV" | trim_cr | while IFS=',' read -r SampleID Index Guide Amplicon; do
  # Basic checks
  if [[ -z "${SampleID:-}" || -z "${Guide:-}" || -z "${Amplicon:-}" ]]; then
    echo "Skipping malformed line: $SampleID"
    continue
  fi

  # Enforce numeric range
  if [[ "$SampleID" =~ ^[0-9]+$ ]]; then
    if (( SampleID < START || SampleID > END )); then
      echo "⏭️  Skip $SampleID (outside ${START}-${END})"
      continue
    fi
  else
    echo "⏭️  Skip non-numeric SampleID '$SampleID' due to active numeric range"
    continue
  fi

  # Input FASTQs expected to be named like: 211_read1.fastq.gz etc.
  R1="${SampleID}_read1.fastq.gz"
  R2="${SampleID}_read2.fastq.gz"

  if [[ ! -s "$R1" || ! -s "$R2" ]]; then
    echo "⚠️  Missing reads for $SampleID (expected $R1 / $R2). Skipping."
    continue
  fi

  # Confirm FASTA exists
  if [[ ! -f "$FASTAS_DIR/${SampleID}.fasta" ]]; then
    echo "⚠️  Missing FASTA for $SampleID. Skipping."
    continue
  fi

  # Temp / intermediate names
  P1="${SampleID}processed.1.fastq.gz"
  P2="${SampleID}processed.2.fastq.gz"
  TRIMMED_R1="${SampleID}_trimmed.1.fastq.gz"
  TRIMMED_R2="${SampleID}_trimmed.2.fastq.gz"
  SAM="${SampleID}.sam"
  BAM="${SampleID}.bam"
  SORTED_BAM="sorted${SampleID}.bam"
  DEDUP_BAM="${SampleID}dedup.bam"
  FINAL_R1="${OUTPUT_FASTQ_DIR}/${SampleID}_R1.fastq.gz"
  FINAL_R2="${OUTPUT_FASTQ_DIR}/${SampleID}_R2.fastq.gz"

  # 1) Extract UMIs
  umi_tools extract \
    -I "$R1" \
    --bc-pattern=NNNNNNNNN \
    --read2-in="$R2" \
    --stdout="$P1" \
    --read2-out="$P2"

  # 2) Adapter trimming
  fastp \
    -i "$P1" \
    -I "$P2" \
    -o "$TRIMMED_R1" \
    -O "$TRIMMED_R2" \
    --adapter_sequence     AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    --adapter_sequence_r2  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --thread "$THREADS" \
    --html "${SampleID}_fastp.html" \
    --json "${SampleID}_fastp.json"

  # 3) Align
  IDX_PREFIX="$FASTAS_DIR/${SampleID}_idx_${SLURM_ARRAY_TASK_ID:-$$}"
  [[ -f "${IDX_PREFIX}.1.bt2" ]] || bowtie2-build -f "$FASTAS_DIR/${SampleID}.fasta" "$IDX_PREFIX"
  bowtie2 -x "$IDX_PREFIX" -1 "$TRIMMED_R1" -2 "$TRIMMED_R2" -S "$SAM" -p "$THREADS"

  # 4) BAM processing
  samtools view -@ "$THREADS" -bS "$SAM" -o "$BAM"
  samtools sort -@ "$THREADS" "$BAM" -o "$SORTED_BAM"
  samtools index "$SORTED_BAM"

  # 5) Deduplication
  umi_tools dedup --paired --method=unique \
    --stdin="$SORTED_BAM" --stdout="$DEDUP_BAM"

  # 6) Convert back to FASTQ
  samtools sort -n -@ "$THREADS" -o "${DEDUP_BAM%.bam}.name.bam" "$DEDUP_BAM"
  samtools fastq -@ "$THREADS" "${DEDUP_BAM%.bam}.name.bam" \
  -1 "$FINAL_R1" -2 "$FINAL_R2" \
  -0 /dev/null -s /dev/null -n

  # 7) Clean up temps
  rm -f "$SAM" "$BAM" "$P1" "$P2" "$TRIMMED_R1" "$TRIMMED_R2"
  rm -f "${IDX_PREFIX}".*.bt2
  echo "✅ Finished $SampleID"
done

echo "🎉 All done!"
