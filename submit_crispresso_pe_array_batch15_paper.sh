#!/bin/bash
#SBATCH --job-name=CRISPRessoRawDataPaperAnalysisPEBatch15Array
#SBATCH --output=logs/crispresso_pe_%A_%a.out
#SBATCH --error=logs/crispresso_pe_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --partition=standard           # change to free if core units are running low
#SBATCH --account=alapinai_lab
#SBATCH --cpus-per-task=4
#SBATCH --array=457-492%10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbking5@asu.edu

mkdir -p logs

# Load environment
module load mamba/latest
source activate ampseq_pipeline

# TSV file (must be in same directory you submit from)
TSV="rawlargebatch15_crispresso_batchfile.tsv"

# Extract this task's sample ID
ROW=$(awk -F'\t' -v id="${SLURM_ARRAY_TASK_ID}" '$2 == id"_R1.fastq.gz"' "$TSV")

# Parse columns (tr -d '\r' strips Windows carriage returns)
NAME=$(echo "$ROW"             | awk -F'\t' '{print $1}'  | tr -d '\r')
FASTQ_R1=$(echo "$ROW"         | awk -F'\t' '{print $2}'  | tr -d '\r')
FASTQ_R2=$(echo "$ROW"         | awk -F'\t' '{print $3}'  | tr -d '\r')
AMPLICON=$(echo "$ROW"         | awk -F'\t' '{print $4}'  | tr -d '\r')
SPACER=$(echo "$ROW"           | awk -F'\t' '{print $5}'  | tr -d '\r')
SCAFFOLD=$(echo "$ROW"         | awk -F'\t' '{print $6}'  | tr -d '\r')
EXTENSION_SEQ=$(echo "$ROW"    | awk -F'\t' '{print $7}'  | tr -d '\r')
RTT_WINDOW=$(echo "$ROW"       | awk -F'\t' '{print $8}'  | tr -d '\r')
PE_OVERRIDE_SEQ=$(echo "$ROW"  | awk -F'\t' '{print $9}'  | tr -d '\r')
NICKING_GUIDE=$(echo "$ROW"    | awk -F'\t' '{print $10}' | tr -d '\r')

# Skip samples with no scaffold — those are base editor samples
if [[ -z "$SCAFFOLD" ]]; then
  echo "⏭️  Skipping $NAME — no scaffold sequence (base editor sample, use BE script instead)"
  exit 0
fi

# Validate required PE parameters are present
if [[ -z "$EXTENSION_SEQ" ]]; then
  echo "❌ ERROR: $NAME — extension_seq is empty in TSV. Please fill in column 7."
  exit 1
fi
if [[ -z "$RTT_WINDOW" ]]; then
  echo "❌ ERROR: $NAME — rtt_window_size is empty in TSV. Please fill in column 8."
  exit 1
fi
if [[ -z "$PE_OVERRIDE_SEQ" ]]; then
  echo "❌ ERROR: $NAME — pe_override_seq is empty in TSV. Please fill in column 9."
  exit 1
fi

echo "Running CRISPResso (Prime Editing) for: $NAME"
echo "  R1: $FASTQ_R1 | R2: $FASTQ_R2"
echo "  Spacer: $SPACER"
echo "  Scaffold: $SCAFFOLD"
echo "  Extension: $EXTENSION_SEQ"
echo "  RTT window: $RTT_WINDOW"
echo "  PE override seq: $PE_OVERRIDE_SEQ"
if [[ -n "$NICKING_GUIDE" ]]; then
  echo "  Nicking guide: $NICKING_GUIDE (PE3 mode)"
else
  echo "  Nicking guide: none (PE2 mode)"
fi

# Build the nicking guide flag conditionally — only added if column 10 is non-empty
NICKING_FLAG=()
if [[ -n "$NICKING_GUIDE" ]]; then
  NICKING_FLAG=(--prime_editing_nicking_guide_seq "$NICKING_GUIDE")
fi

CRISPResso \
  --fastq_r1 "$FASTQ_R1" \
  --fastq_r2 "$FASTQ_R2" \
  --amplicon_seq "$AMPLICON" \
  --prime_editing_pegRNA_spacer_seq "$SPACER" \
  --prime_editing_pegRNA_extension_seq "$EXTENSION_SEQ" \
  --prime_editing_pegRNA_scaffold_seq "$SCAFFOLD" \
  --name "$NAME" \
  --plot_window_size 40 \
  --prime_editing_pegRNA_extension_quantification_window_size "$RTT_WINDOW" \
  --prime_editing_pegRNA_scaffold_min_match_length 1 \
  --min_bp_quality_or_N 0 \
  --exclude_bp_from_left 0 \
  --exclude_bp_from_right 0 \
  --assign_ambiguous_alignments_to_first_reference \
  --prime_editing_override_prime_edited_ref_seq "$PE_OVERRIDE_SEQ" \
  "${NICKING_FLAG[@]}" \
  --output_folder "CRISPRessoBatch15PaperAnalysis_PEsamples"

echo "✅ Finished $NAME"
