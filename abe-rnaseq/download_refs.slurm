#!/bin/bash
#SBATCH -J dl_refs
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=madelk1@hs.uci.edu
#SBATCH -A madelk1
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH -o /pub/madelk1/rnaseq/logs/%x-%j.out
#SBATCH -e /pub/madelk1/rnaseq/logs/%x-%j.err

set -euo pipefail
module purge
module load samtools/1.15.1 java/17 picard-tools/3.3.0

REF=/pub/madelk1/rnaseq/refs/hg38
mkdir -p "$REF" && cd "$REF"

# Genome + GTF
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip -f GRCh38.primary_assembly.genome.fa.gz
gunzip -f gencode.v44.annotation.gtf.gz

# FASTA indices
samtools faidx GRCh38.primary_assembly.genome.fa
JAVA=$(which java)
PICARD_JAR=/opt/apps/picard-tools/3.3.0/picard.jar
"$JAVA" -Xmx4g -jar "$PICARD_JAR" CreateSequenceDictionary \
  R=GRCh38.primary_assembly.genome.fa \
  O=GRCh38.primary_assembly.genome.dict

# dbSNP (check periodically for updates — this file gets revised)
wget -c -q --show-progress \
  https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -c -q --show-progress \
  https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
