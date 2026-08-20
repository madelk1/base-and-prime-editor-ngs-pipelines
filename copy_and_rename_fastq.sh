#!/bin/bash

# Create output directory
mkdir -p "Renamed_Files"

# Loop over all matching FASTQ files like: 75_R1_001.fastq.gz
for file in *_R[12]_001.fastq.gz; do
    # Sample ID = text before first underscore (e.g. "75")
    sample="${file%%_*}"

    # Determine read type
    if [[ "$file" == *_R1_001.fastq.gz ]]; then
        read_type="read1"
    elif [[ "$file" == *_R2_001.fastq.gz ]]; then
        read_type="read2"
    else
        echo "Unknown read type in file: $file"
        continue
    fi

    # New filename: e.g. 75_read1.fastq.gz
    new_name="${sample}_${read_type}.fastq.gz"

    # Copy (or mv) to Renamed_Files
    cp "$file" "Renamed_Files/$new_name"
    echo "Copied $file -> Renamed_Files/$new_name"
done
