#!/bin/bash

# Create a logs folder if it doesn't already exist
mkdir -p ./logs

GENOME=/projects/micb405/resources/genomes/mouse/mm10/bwa_index/mm10.fa
CHIPDATA=/projects/micb405/data/mouse/project_1/ChIP-seq/

declare -a arr=("Naive_H3K27ac" "Naive_Input" "Primary_H3K27ac" "Primary_Input")

# Align the fastq files
for i in "${arr[@]}"
do
    echo "Starting bwa for $i experimental"
    bwa mem -t 16 $GENOME $CHIPDATA/${i}_1.fastq $CHIPDATA/${i}_2.fastq >./${i}.sam 2>logs/bwa_${i}.log
done
