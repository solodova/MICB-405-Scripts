#!/bin/bash

ZIPPED_GENOME=/projects/micb405/resources/project_2/2019/SaanichInlet_200m/MetaT_reads/7753.4.82794.ATTCCT.qtrim.3ptrim.artifact.rRNA.clean.fastq.gz
SCRIPTS_DIR=/projects/micb405/project1/Team7/Proj2_scripts

echo "Uncompressing:
 ${ZIPPED_GENOME}"

gunzip -c ${ZIPPED_GENOME} > ${SCRIPTS_DIR}/MetaT_reads.fastq

echo "Uncompressing complete!"
