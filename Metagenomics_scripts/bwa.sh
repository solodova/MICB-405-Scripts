#!/bin/bash

GENOME=/projects/micb405/project1/Team7/Proj2_scripts/MetaT_reads.fastq
PROKKA_DIR=/projects/micb405/project1/Team7/Prokka_output

dir1=${PROKKA_DIR}/SaanichInlet_200m_${1}

if [[ ! -d "${dir1}" ]]; then
  echo "First argument is invalid! ${dir1} MAG not found."
  exit 1
fi

mkdir -p logs

ref1Count=`ls ${dir1} -1 | grep .ffn | wc -l`

ref1=${dir1}/`ls ${dir1} -1 | grep .ffn\$`

if [[ ${ref1Count} -eq 1 ]]; then
    echo "Indexing the following file:
    ${ref1}"
    bwa index ${ref1}
else
    echo "${1} is already indexed! Moving on..."
fi

echo "Aligning the reference files to the read"
echo "ref1 = ${ref1}"
echo "GENOME = ${GENOME}"

bwa mem -t ${2:-16} ${ref1} ${GENOME} > ./${1}.sam 2> logs/bwa_${1}.log &
