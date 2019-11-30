#!/bin/bash
MAG_PATH=/projects/micb405/resources/project_2/2019/SaanichInlet_200m/MetaBAT2_SaanichInlet_200m/MedQPlus_MAGs/

ARCHAEA=("113" "137" "148" "154" "174" "189" "201" "210" "220" "53" "55" "56" "57" "63" "65" "70")

for MAG in ${MAG_PATH}* ; do
	MAG_NAME=$(basename ${MAG%.*})
	MAG_NUM=${MAG_NAME##*.}
	KINGDOM=Bacteria
	if [[ " ${ARCHAEA[@]} " =~ " ${MAG_NUM} " ]]; then
		KINGDOM=Archaea
	fi
	PROKKA_OUTDIR="Prokka_output/SaanichInlet_200m_${MAG_NUM}_${KINGDOM}/"
	echo "Running prokka, outputting to" $PROKKA_OUTDIR
	prokka \
	--outdir $PROKKA_OUTDIR \
	--prefix "MAG${MAG_NUM}" \
	--kingdom $KINGDOM \
	--cpus 2 \
	$MAG
	wait ;
done ;
