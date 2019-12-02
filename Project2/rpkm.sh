#!/bin/bash

# path to RPKM executable
RPKM=/projects/micb405/resources/project_2/2019/rpkm
# folder containing Prokka output for MAGs
PROKKA_DIR=/projects/micb405/project1/Team7/Prokka_output

# this script takes in a command line argument specifying the ID of the MAG,
# and uses this ID to grab MAG-specific SAM files (from bwa mem) and .ffn files (from Prokka)
# and finally uses these to run rpkm.


dir1=${PROKKA_DIR}/SaanichInlet_200m_${1}

if [[ ! -d "${dir1}" ]]; then
  echo "First argument is invalid! ${dir1} MAG not found."
  exit 1
fi

ffn1=${dir1}/`ls ${dir1} -1 | grep .ffn\$`

echo "executing RPKM for ${1}"

${RPKM} \
-c ${ffn1} \
-a ./${1}.sam \
-o ./${1}.csv \
--verbose &
