#!/bin/bash

declare -a arr=("Naive" "Primary")

for i in "${arr[@]}"
do
    echo "running macs2 for $i cells"
    macs2 callpeak \
    -t ${i}_H3K27ac.sorted.mkdup.bam \
    -c ${i}_Input.sorted.mkdup.bam \
    -f BAMPE -g mm -n ${i}_H3K27ac -B -q 0.01 \
    2> logs/${i}_callpeak.log
done
