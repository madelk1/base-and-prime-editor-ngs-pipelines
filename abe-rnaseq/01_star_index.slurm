#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=madelk1@hs.uci.edu
#SBATCH -A madelk1
#SBATCH -J star_index
#SBATCH -p standard
#SBATCH -t 04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -o /pub/madelk1/rnaseq/logs/star_index-%j.out
#SBATCH -e /pub/madelk1/rnaseq/logs/star_index-%j.err

module purge
module load star/2.7.10a

IDX=/pub/madelk1/rnaseq/refs/hg38/star_index
FA=/pub/madelk1/rnaseq/refs/hg38/GRCh38.primary_assembly.genome.fa
GTF=/pub/madelk1/rnaseq/refs/hg38/gencode.v44.annotation.gtf
mkdir -p "$IDX"

STAR --runThreadN ${SLURM_CPUS_PER_TASK} --runMode genomeGenerate \
     --genomeDir "$IDX" \
     --genomeFastaFiles "$FA" \
     --sjdbGTFfile "$GTF" \
     --sjdbOverhang 149
