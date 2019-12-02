#!/bin/bash

SUMMITS=_H3K27ac_summits
PEAKS=_H3K27ac_peaks

# Need to transform bed files and narrowPeak files to be compatible
# with Bedtools. Summit files need 4 columns and narrow peaks need 6

bedCount=`ls -1 | grep .4col.bed | wc -l`

declare -a arr=("Naive" "Primary")

if [[ ${bedCount} == 0 ]]
    then
        echo "Transforming summit files"
        for i in "${arr[@]}"
        do
        	cut -f1-4 ${i}${SUMMITS}.bed > ${i}${SUMMITS}.4col.bed
        done
fi

peakCount=`ls -1 | grep .6col.bed | wc -l`

if [[ ${peakCount} == 0 ]]
    then
        echo "Transforming narrowPeak files"
        for i in "${arr[@]}"
        do
        	cut -f1-6 ${i}${PEAKS}.narrowPeak > ${i}${PEAKS}.6col.bed
        done
fi

# Use bedtools subtract to difference between the two types

echo "Running bedtools for summit files"

bedtools subtract \
    -a ${arr[1]}${SUMMITS}.4col.bed \
    -b ${arr[0]}${SUMMITS}.4col.bed \
    > summit_difference.txt

echo "Running bedtools for narrowPeak files"

bedtools subtract \
    -a ${arr[1]}${PEAKS}.6col.bed \
    -b ${arr[0]}${PEAKS}.6col.bed \
    > narrowPeak_difference.txt
