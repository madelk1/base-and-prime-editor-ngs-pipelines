#!/bin/bash
#SBATCH --job-name=haplotypecaller
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=alapinai_lab
#SBATCH -o /pub/madelk1/HiFiRNASeq/logs/haplotypecaller-%j.out
#SBATCH -e /pub/madelk1/HiFiRNASeq/logs/haplotypecaller-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madelk1@uci.edu

set -e
set -o pipefail

sample=$1

# Paths
fasta="/pub/madelk1/rnaseq/refs/hg38/GRCh38.primary_assembly.genome.fa"
bam_dir="/pub/madelk1/HiFiRNASeq/forVariantCalling"
output_dir="/pub/madelk1/HiFiRNASeq/variantCalls"

mkdir -p ${output_dir}

# Load GATK4
module load gatk/4.6.2.0

echo "========================================="
echo "Processing sample: ${sample}"
echo "Started at: $(date)"
echo "========================================="

# Run HaplotypeCaller
gatk --java-options "-Xmx60G" HaplotypeCaller \
 -R ${fasta} \
 -I "${bam_dir}/${sample}.bam" \
 -O "${output_dir}/${sample}.vcf.gz" \
 --dont-use-soft-clipped-bases \
 --standard-min-confidence-threshold-for-calling 20.0 \
 --native-pair-hmm-threads 16

echo "Tool returned: $?"

echo "========================================="
echo "Finished at: $(date)"
echo "Output VCF: ${output_dir}/${sample}.vcf.gz"
echo "========================================="
