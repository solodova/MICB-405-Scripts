#!/bin/bash
RNA_PREFIX=/projects/micb405/data/mouse/project_1/RNA-seq
# #### O: Path to group 7 folder
GROUP7=/projects/micb405/project1/Team7
#mkdir $GROUP7/FastQC_output/
#
# #### O: all the RNA SEQ sets (read 1) for primary/activated T cells
RNA_SEQ_SET_NAIVE=($RNA_PREFIX/Naive_1.1.fastq $RNA_PREFIX/Naive_1.2.fastq $RNA_prefix/Naive_1.3.fastq)
RNA_SEQ_SET_PRIMARY=($RNA_PREFIX/Primary_1.1.fastq $RNA_PREFIX/Primary_1.2.fastq $RNA_PREFIX/Primary_1.3.fastq)

for fastq in ${RNA_SEQ_SET_PRIMARY[@]}; do
  PREFIX=$(basename $fastq .fastq) ;
  CURR_ID=${PREFIX#Primary_1.} ;
  CURR_FASTQ=Primary_?.$CURR_ID.fastq ;
  echo $PREFIX $CURR_ID $CURR_FASTQ ;
  echo "Doing FastQC for " $RNA_PREFIX/$CURR_FASTQ ;
  fastqc --threads 1 \
  -o $GROUP7/FastQC_output/PRIMARY$CURR_ID/ \
  $RNA_PREFIX/$CURR_FASTQ ;
  wait
done

