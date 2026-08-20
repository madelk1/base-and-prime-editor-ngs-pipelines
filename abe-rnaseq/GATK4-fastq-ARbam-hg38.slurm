#!/bin/bash
#SBATCH --job-name=gatk4_pipeline_hifi
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --partition=standard
#SBATCH --account=alapinai_lab
#SBATCH -o /pub/madelk1/HiFiRNASeq/logs/gatk4_pipeline-%j.out
#SBATCH -e /pub/madelk1/HiFiRNASeq/logs/gatk4_pipeline-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madelk1@uci.edu
set -e
set -o pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch $0 <sample_name>"
  echo "Example: sbatch $0 001_CTRL1"
  exit 1
fi

sample=$1

BASEDIR="/pub/madelk1/HiFiRNASeq"
REF="/pub/madelk1/rnaseq/refs/hg38/GRCh38.primary_assembly.genome.fa"
STARINDEX="/pub/madelk1/rnaseq/refs/hg38/star_index"
GTF="/pub/madelk1/rnaseq/refs/hg38/gencode.v44.annotation.gtf"
DBSNP="/pub/madelk1/rnaseq/refs/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"

mkdir -p ${BASEDIR}/aligned
mkdir -p ${BASEDIR}/temp
mkdir -p ${BASEDIR}/metrics
mkdir -p ${BASEDIR}/forVariantCalling
mkdir -p ${BASEDIR}/logs

cd ${BASEDIR}

module load star/2.7.10a
module load gatk/4.6.2.0
module load samtools/1.15.1

#==========================================
# Step 1: STAR alignment
# ADDED: --outSAMmapqUnique 60 (GATK4 best practice)
#==========================================
echo "[$(date)] STAR alignment..."
STAR --runThreadN 16 \
 --genomeDir ${STARINDEX} \
 --sjdbGTFfile ${GTF} \
 --readFilesIn ${sample}.R1.fastq.gz ${sample}.R2.fastq.gz \
 --readFilesCommand zcat \
 --outFileNamePrefix aligned/${sample} \
 --outSAMtype BAM SortedByCoordinate \
 --outFilterMultimapNmax 1 \
 --outSAMmapqUnique 60 \
 --twopassMode Basic \
 --quantMode GeneCounts
samtools index aligned/${sample}Aligned.sortedByCoord.out.bam

#==========================================
# Step 2: AddOrReplaceReadGroups
# FIXED: RGSM ${sample} (not "1")
#==========================================
echo "[$(date)] Adding read groups..."
gatk AddOrReplaceReadGroups \
 -I aligned/${sample}Aligned.sortedByCoord.out.bam \
 -O temp/${sample}.temp1.bam \
 -SO coordinate \
 -RGID ${sample} \
 -RGLB lib1 \
 -RGPL illumina \
 -RGPU unit1 \
 -RGSM ${sample}
echo "Tool returned: $?"

#==========================================
# Step 3: MarkDuplicates
# KEPT: REMOVE_DUPLICATES true (correct for RNA-seq)
#==========================================
echo "[$(date)] Marking duplicates..."
gatk MarkDuplicates \
 -I temp/${sample}.temp1.bam \
 -O temp/${sample}.temp2.bam \
 -M metrics/${sample}.rmdup.log \
 --REMOVE_DUPLICATES true \
 --ASSUME_SORTED true \
 --VALIDATION_STRINGENCY LENIENT
echo "Tool returned: $?"
samtools index temp/${sample}.temp2.bam

#==========================================
# Step 4: SplitNCigarReads
# REMOVED: -rf ReassignOneMappingQuality (handled by STAR now)
#==========================================
echo "[$(date)] SplitNCigarReads..."
gatk SplitNCigarReads \
 -R ${REF} \
 -I temp/${sample}.temp2.bam \
 -O temp/${sample}.temp3.bam
echo "Tool returned: $?"

#==========================================
# Step 5: BaseRecalibrator - Pass 1
#==========================================
echo "[$(date)] BaseRecalibrator pass 1..."
gatk BaseRecalibrator \
 -R ${REF} \
 -I temp/${sample}.temp3.bam \
 -O metrics/${sample}.recal_data.table \
 --known-sites ${DBSNP}
echo "Tool returned: $?"

#==========================================
# Step 6: ApplyBQSR
# FIXED: --bqsr-recal-file (full parameter name)
#==========================================
echo "[$(date)] ApplyBQSR..."
gatk ApplyBQSR \
 -R ${REF} \
 -I temp/${sample}.temp3.bam \
 --bqsr-recal-file metrics/${sample}.recal_data.table \
 -O forVariantCalling/${sample}.bam
echo "Tool returned: $?"

#==========================================
# Step 7: BaseRecalibrator - Pass 2 (QC)
# FIXED: Run on recalibrated BAM (not temp3.bam)
#==========================================
echo "[$(date)] BaseRecalibrator pass 2..."
gatk BaseRecalibrator \
 -R ${REF} \
 -I forVariantCalling/${sample}.bam \
 -O metrics/${sample}.post_recal_data.table \
 --known-sites ${DBSNP}
echo "Tool returned: $?"

# Index final BAM
samtools index forVariantCalling/${sample}.bam

echo "========================================="
echo "Pipeline complete for ${sample}!"
echo "Finished at: $(date)"
echo "========================================="
echo ""
echo "Final BAM: ${BASEDIR}/forVariantCalling/${sample}.bam"

chmod +x /pub/madelk1/HiFiRNASeq/GATK4updatedscripts/GATK4-fastq-ARbam-hg38.slurm
