#!/bin/bash

# Get relative depth of fragments for naive and primary cells

declare -a arr=("Naive" "Primary")
FRAG_STRING="fragments after filtering in control"

for i in "${arr[@]}"
do
    eval "Depth_${i}=$(grep "$FRAG_STRING" ./logs/${i}_callpeak.log | awk '{print $NF}')"
done

TREATMENT="_H3K27ac_treat_pileup"
CONTROL="_H3K27ac_control_lambda"
PREFIX="diff_${arr[1]}_vs_${arr[0]}"

# Find differentially expressed genes between naive and primary cells

macs2 bdgdiff \
--t1 ${arr[1]}${TREATMENT}.bdg \
--c1 ${arr[1]}${CONTROL}.bdg \
--t2 ${arr[0]}${TREATMENT}.bdg \
--c2 ${arr[0]}${CONTROL}.bdg \
--d1 ${Depth_Primary} \
--d2 ${Depth_Naive} -g 60 -l 120 \
--o-prefix ${PREFIX}

# Format the output files so they don't have headers

declare -a arr2=("cond1" "cond2" "common")

for i in "${arr2[@]}"
do
    tail -n +2 ${PREFIX}_c3.0_${i}.bed > ${PREFIX}_${i}_formatted.bed
done
