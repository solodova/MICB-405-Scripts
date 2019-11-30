#!/bin/bash

# Create output dir for fastqc if it doesn't already exist
mkdir -p ./FastQC_output

CHIPDATA=/projects/micb405/data/mouse/project_1/ChIP-seq/

declare -a arr=("Naive_H3K27ac" "Naive_Input" "Primary_H3K27ac" "Primary_Input")

# Align the fastq files
for i in "${arr[@]}"
do
    echo "Starting bwa for $i experimental"
    fastqc --threads 2 -o ./FastQC_output/ ${CHIPDATA}/${i}_1.fastq
done
