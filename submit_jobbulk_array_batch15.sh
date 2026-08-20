#!/bin/bash
#SBATCH --job-name=array_pipeline_batch15
#SBATCH --output=logs/pipeline_%A_%a.out
#SBATCH --error=logs/pipeline_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --partition=standard #change to free if core units are running low
#SBATCH --account=alapinai_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=1-456,493-519%20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mbking5@asu.edu

# Ensure logs directory exists
mkdir -p logs

# Load necessary modules
module load mamba/latest

# Pass the array task ID as a single-sample range (e.g. "5-5") to the original pipeline script
bash batch_15_pipeline_script.sh "${SLURM_ARRAY_TASK_ID}-${SLURM_ARRAY_TASK_ID}"
