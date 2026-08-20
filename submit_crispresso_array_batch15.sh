#!/bin/bash
#SBATCH --job-name=CRISPRessoBatch14Array
#SBATCH --output=logs/crispresso_%A_%a.out
#SBATCH --error=logs/crispresso_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --partition=standard #change to free if core units are running low
#SBATCH --account=alapinai_lab
#SBATCH --cpus-per-task=2
#SBATCH --array=1-456,493-519%20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbking5@asu.edu

mkdir -p logs

# Load environment
module load mamba/latest
source activate ampseq_pipeline

# TSV file (must be in the same directory you submit from)
TSV="processedlargebatch15_crispresso_batchfile.tsv"

# Extract this task's sampleID (skip header, grab line matching task ID)
ROW=$(awk -F'\t' -v id="${SLURM_ARRAY_TASK_ID}" '$2 == id"_R1.fastq.gz"' "$TSV")

# Parse columns (tr -d '\r' strips Windows carriage returns)
NAME=$(echo "$ROW"     | awk -F'\t' '{print $1}' | tr -d '\r')
FASTQ_R1=$(echo "$ROW" | awk -F'\t' '{print $2}' | tr -d '\r')
FASTQ_R2=$(echo "$ROW" | awk -F'\t' '{print $3}' | tr -d '\r')
AMPLICON=$(echo "$ROW" | awk -F'\t' '{print $4}' | tr -d '\r')
GUIDE=$(echo "$ROW"    | awk -F'\t' '{print $5}' | tr -d '\r')


echo "Running CRISPResso for: $NAME"
echo "  R1: $FASTQ_R1 | R2: $FASTQ_R2"
echo "  Guide: $GUIDE"

CRISPResso \
  --fastq_r1 "$FASTQ_R1" \
  --fastq_r2 "$FASTQ_R2" \
  --amplicon_seq "$AMPLICON" \
  --guide_seq "$GUIDE" \
  --name "$NAME" \
  --quantification_window_size 20 \
  --quantification_window_center -10 \
  --exclude_bp_from_left 0 \
  --exclude_bp_from_right 0 \
  --base_editor_output \
  --output_folder "CRISPRessoBatch15Array"

echo "✅ Finished $NAME"
