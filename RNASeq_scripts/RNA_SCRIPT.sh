#!/bin/bash
# # First make a directory to download and store the reference files,
# # dont worry about downloading the assembly.fa but download the gtf file
#
# #mkdir $HOME/References/STAR_tutorial_mus_musculus_GRCm38/
# #cd $HOME/References/STAR_tutorial_mus_musculus_GRCm38/
#
# #### O: we already have the reference genome and annotation file
# #wget ftp://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
# #wget ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz
#
# # Now we can create an index on the assemby
# # Note the file extension is gone, this is because unzipped the files
#
# #### O: paths to genome index files, fasta file, annotation file respectively
# GENOME_DIR=/projects/micb405/resources/STAR_tutorial/STAR_index_musmusculus_mm10/
# GENOME_FASTA=$HOME/References/STAR_tutorial_mus_musculus_GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa
GENOME_GTF=/projects/micb405/resources/STAR_tutorial/STAR_index_musmusculus_mm10/mm10.gtf
#
# #### O: we don't need to run this because the index was generated for us
# # STAR --runThreadN 30 \
# # --runMode genomeGenerate \
# # --genomeDir $GENOME_DIR \
# # --genomeFastaFiles $GENOME_FASTA \
# # --sjdbGTFfile $GENOME_GTF \
# # --sjdbOverhang 74 & \
#
# #### O: RNA SEQ datasets are paired-end
# #### O: fastq file naming convention:  <type_read.ID/sample>.fastq
#
# #### O: Path to RNA seq data folder
# RNA_PREFIX=/projects/micb405/data/mouse/project_1/RNA-seq/
# #### O: Path to group 7 folder
GROUP7=/projects/micb405/project1/Team7/
STAR_OUTPUT=$GROUP7/STAR_alignment/
#
# #### O: all the RNA SEQ sets (read 1) for primary/activated T cells
# RNA_SEQ_SET_PRIMARY=($RNA_PREFIX/Primary_1.1.fastq $RNA_PREFIX/Primary_1.2.fastq $RNA_PREFIX/Primary_1.3.fastq)
# # O: iterate over RNA SEQ sets
# for i in ${RNA_SEQ_SET_PRIMARY[@]}; do
#   # basename will get file name from path
#   # eg. OUTPREFIX for first RNASEQ set is Primary_1.1
#   OUTPREFIX=$(basename $i .fastq ) ;
#   # set the prefix for the output file to be PRIMARY<ID>
#   OUTPREFIX=PRIMARY${OUTPREFIX#Primary_1.} ;
#   echo "Starting alignment for" $OUTPREFIX ;
#   STAR --genomeDir $GENOME_DIR \
#   --runThreadN 1 \
#   --readFilesIn $i ${i/Primary_1/Primary_2} \
#   --outFileNamePrefix $STAR_OUTPUT/$OUTPREFIX \
#   --outSAMtype BAM SortedByCoordinate \
#   --outSAMunmapped Within \
#   --outSAMattributes Standard ;
#   wait
# done
#
# #all the RNA SEQ sets (read 1) for naive T cells
# RNA_SEQ_SET_NAIVE=($RNA_PREFIX/Naive_1.1.fastq $RNA_PREFIX/Naive_1.2.fastq $RNA_PREFIX/Naive_1.3.fastq)
# for i in ${RNA_SEQ_SET_NAIVE[@]}; do
#   # basename will get file name from path
#   # eg. OUTPREFIX for first RNASEQ set is Naive_1.1
#   OUTPREFIX=$(basename $i .fastq ) ;
#   # set the prefix for the output file to be NAIVE<ID>
#   OUTPREFIX=NAIVE${OUTPREFIX#Naive_1.} ;
#   echo "Starting alignment for" $OUTPREFIX ;
#   STAR --genomeDir $GENOME_DIR \
#   --runThreadN 1 \
#   --readFilesIn $i ${i/Naive_1/Naive_2} \
#   --outFileNamePrefix $STAR_OUTPUT/$OUTPREFIX \
#   --outSAMtype BAM SortedByCoordinate \
#   --outSAMunmapped Within \
#   --outSAMattributes Standard ;
#   wait
# done

#if your files are gzipped you can use the command below to run zcat to unzip the files
#--readFilesCommand zcat \
#Once this is done lets check out the alignment results in our alignment folder (Log.final.out)

###Do quality control

RNA_BAM_P1=$STAR_OUTPUT/PRIMARY1Aligned.sortedByCoord.out.bam
RNA_BAM_P2=$STAR_OUTPUT/PRIMARY2Aligned.sortedByCoord.out.bam
RNA_BAM_P3=$STAR_OUTPUT/PRIMARY3Aligned.sortedByCoord.out.bam
RNA_BAM_N1=$STAR_OUTPUT/NAIVE1Aligned.sortedByCoord.out.bam
RNA_BAM_N2=$STAR_OUTPUT/NAIVE2Aligned.sortedByCoord.out.bam
RNA_BAM_N3=$STAR_OUTPUT/NAIVE3Aligned.sortedByCoord.out.bam

RNA_BAMS=($RNA_BAM_P1 $RNA_BAM_P2 $RNA_BAM_P3 $RNA_BAM_N1 $RNA_BAM_N2 $RNA_BAM_N3)

#for bam in ${RNA_BAMS[@]}; do
#  OUTPREFIX=$(basename $bam Aligned.sortedByCoord.out.bam)
#  OUTNAME=htseq_counts_$OUTPREFIX.txt
#  htseq-count \
#  -f bam \
#  -m union \
#  -i gene_id \
#  -r pos \
#  -s reverse \
#  $bam  \
#  $GENOME_GTF > $OUTNAME ;
#  wait
#done

# ###get counts using htseq
# htseq-count \
# -f bam \
# -m union \
# -i gene_id \
# -r pos \
# -s reverse \
# ${RNA_BAMS[@]}  \
# $GENOME_GTF > htseq_counts_star.txt
# wait

#cd ~ ; wget http://hartleys.github.io/QoRTs/QoRTs-STABLE.jar ;

for bam in ${RNA_BAMS[@]}; do
  OUTPREFIX=$(basename $bam Aligned.sortedByCoord.out.bam)
  QC_DIR=$GROUP7/STAR_qc/$OUTPREFIX
  mkdir $QC_DIR
  echo "Starting QC for" $QC_DIR ;
  java -jar $HOME/QoRTs-STABLE.jar QC \
  --generatePlots \
  --stranded \
  $bam \
  $GENOME_GTF \
  $QC_DIR ;
  wait
done
