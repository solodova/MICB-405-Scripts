#!/bin/bash

declare -a arr=("Naive_H3K27ac" "Naive_Input" "Primary_H3K27ac" "Primary_Input")

for i in "${arr[@]}"
do
    echo "Running sambamba for $i"
    sambamba view -S -f bam -o $i.bam $i.sam                # Convert to bam format
    sambamba sort -t 8 $i.bam                               # Sort bam files
    sambamba markdup -t 8 $i.sorted.bam $i.sorted.mkdup.bam # Mark duplicates
    echo "$i is ready for macs2, deleting intermediate files"
    rm $i.bam $i.sam $i.sorted.bam
done

echo "Removing .bai files"
rm *.bai
