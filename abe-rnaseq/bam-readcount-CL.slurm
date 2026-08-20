#!/bin/bash
#SBATCH --job-name=bam_readcount
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=72:00:00
#SBATCH --partition=standard
#SBATCH --account=alapinai_lab
#SBATCH -o /pub/madelk1/HiFiRNASeq/logs/bam_readcount-%j.out
#SBATCH -e /pub/madelk1/HiFiRNASeq/logs/bam_readcount-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madelk1@uci.edu

sample=$1

# Paths
fasta="/pub/madelk1/rnaseq/refs/hg38/GRCh38.primary_assembly.genome.fa"
bam_dir="/pub/madelk1/HiFiRNASeq/forVariantCalling"
output_dir="/pub/madelk1/HiFiRNASeq/bam-readcount"
mkdir -p ${output_dir}

echo "========================================="
echo "Processing sample: ${sample}"
echo "Started at: $(date)"
echo "========================================="

# Activate conda environment
source /opt/apps/anaconda/2020.07/etc/profile.d/conda.sh
conda activate rna_aux

# Run bam-readcount with quality filters
# -w 1 = site filtering (1 = no filtering)
bam-readcount -w 1 -f ${fasta} "${bam_dir}/${sample}.bam" | \
 sed 's/:/\t/g' | \
 awk 'BEGIN{print "CHR\tBP\tREF\tTOTAL\tA\tC\tG\tT\tN"}; {print $1,$2,$3,$4,$20,$34,$48,$62,$76}' OFS='\t' | \
 gzip > "${output_dir}/${sample}.all_positions.txt.gz"

conda deactivate

echo "========================================="
echo "Finished at: $(date)"
echo "Output: ${output_dir}/${sample}.all_positions.txt.gz"
echo "========================================="
