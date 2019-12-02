#!/bin/bash

# path to un-compressed metatranscriptomic reads
GENOME=/projects/micb405/project1/Team7/Proj2_scripts/MetaT_reads.fastq
# path to Prokka output (.ffn files with predicted ORFs)
PROKKA_DIR=/projects/micb405/project1/Team7/Prokka_output

# this script takes in a command line argument specifying ID of MAG to run bwa mem with
# the ID provided is used to specify a MAG-specific .ffn file from Prokka to be used to index the genome,
# which is then used as a reference to align the metatranscriptomic reads

dir1=${PROKKA_DIR}/SaanichInlet_200m_${1}

if [[ ! -d "${dir1}" ]]; then
  echo "First argument is invalid! ${dir1} MAG not found."
  exit 1
fi

mkdir -p logs

# make an index file for specified MAG
ref1Count=`ls ${dir1} -1 | grep .ffn | wc -l`
ref1=${dir1}/`ls ${dir1} -1 | grep .ffn\$`
if [[ ${ref1Count} -eq 1 ]]; then
    echo "Indexing the following file:
    ${ref1}"
    bwa index ${ref1}
else
    echo "${1} is already indexed! Moving on..."
fi

# perform alignment using bwa mem
echo "Aligning the reads to the reference file"
echo "ref1 = ${ref1}"
echo "GENOME = ${GENOME}"
bwa mem -t ${2:-16} ${ref1} ${GENOME} > ./${1}.sam 2> logs/bwa_${1}.log &
