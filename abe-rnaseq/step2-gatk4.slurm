#!/bin/bash
#SBATCH --job-name=step2_merge
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=standard
#SBATCH --account=alapinai_lab
#SBATCH -o /pub/madelk1/HiFiRNASeq/logs/step2_merge-%j.out
#SBATCH -e /pub/madelk1/HiFiRNASeq/logs/step2_merge-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=madelk1@uci.edu

set -e
set -o pipefail

if [[ $# -lt 4 ]]; then
  echo "Usage: sbatch $(basename $0)
  --vcf_control=control_vcf_file.vcf.gz
  --readcounts_control=control_bam_readcounts.txt.gz
  --vcf_treated=treated_vcf_file.vcf.gz
  --readcounts_treated=treated_bam_readcounts.txt.gz
  " >&2
  exit 1
fi

for i in "$@"
do
case $i in
  --vcf_control=*)
    VCF_CONTROL="${i#*=}"
    ;;
  --readcounts_control=*)
    READCOUNTS_CONTROL_GZIP="${i#*=}"
    ;;
  --vcf_treated=*)
    VCF_TREATED="${i#*=}"
    ;;
  --readcounts_treated=*)
    READCOUNTS_TREATED_GZIP="${i#*=}"
    ;;
  -h)
    echo "Usage: sbatch $(basename $0)
  --vcf_control=control_vcf_file.vcf.gz
  --readcounts_control=control_bam_readcounts.txt.gz
  --vcf_treated=treated_vcf_file.vcf.gz
  --readcounts_treated=treated_bam_readcounts.txt.gz
  " >&2
    exit 1
    ;;
  *)
    ;;
esac
done

echo "==========================================="
echo "Step 2: Merging VCF with bam-readcount data"
echo "Started at: $(date)"
echo "==========================================="

READCOUNTS_CONTROL=$(echo ${READCOUNTS_CONTROL_GZIP} | sed 's/\.gz//g')
READCOUNTS_TREATED=$(echo ${READCOUNTS_TREATED_GZIP} | sed 's/\.gz//g')

# Sort files to be joined for control
LC_ALL=C zcat ${READCOUNTS_CONTROL_GZIP} | awk '{print $1"_"$2"\t"$0}' | sort -k1,1 > ${READCOUNTS_CONTROL}.sorted
LC_ALL=C zgrep -v "^#" ${VCF_CONTROL} | awk '{print $1"_"$2"\t"$0}' | sort -k1,1 > ${VCF_CONTROL}.sorted
join ${VCF_CONTROL}.sorted ${READCOUNTS_CONTROL}.sorted -1 1 -2 1 | sed 's/ /\t/g' | cut -f1-7,14-21 > ${VCF_CONTROL}.sorted.vcf_with_coverage.txt

# Sort files to be joined for treated
LC_ALL=C zcat ${READCOUNTS_TREATED_GZIP} | awk '{print $1"_"$2"\t"$0}' | sort -k1,1 > ${READCOUNTS_TREATED}.sorted
LC_ALL=C zgrep -v "^#" ${VCF_TREATED} | awk '{print $1"_"$2"\t"$0}' | sort -k1,1 > ${VCF_TREATED}.sorted
join ${VCF_TREATED}.sorted ${READCOUNTS_TREATED}.sorted -1 1 -2 1 | sed 's/ /\t/g' | cut -f1-7,14-21 > ${VCF_TREATED}.sorted.vcf_with_coverage.txt

# Join treated vcf with control readcounts
join ${VCF_TREATED}.sorted.vcf_with_coverage.txt ${READCOUNTS_CONTROL}.sorted -1 1 -2 1 -a 1 | sed 's/ /\t/g' | awk '{ if (NF == 14) {print $0"\t0\t0\t0\t0\t0\t0\t0\t0\t0"} else {print $0}}' > ${VCF_TREATED}.sorted.vcf_with_coverage.be_and_control.txt

echo "==========================================="
echo "Finished at: $(date)"
echo "Output: ${VCF_TREATED}.sorted.vcf_with_coverage.be_and_control.txt"
echo "==========================================="
