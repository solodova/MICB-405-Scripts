#!/bin/bash
# path to our medium and high quality MAGs
MAG_PATH=/projects/micb405/resources/project_2/2019/SaanichInlet_200m/MetaBAT2_SaanichInlet_200m/MedQPlus_MAGs/

# path to group project folder
GROUP_PATH=/projects/micb405/project1/Team7

# MAG IDs which gtbdk predicted to be Archaea
ARCHAEA=("113" "137" "148" "154" "174" "189" "201" "210" "220" "53" "55" "56" "57" "63" "65" "70")

# Go through all MAGs, run Prokka to get predicted ORFs
for MAG in ${MAG_PATH}* ; do
	MAG_NAME=$(basename ${MAG%.*})
	MAG_NUM=${MAG_NAME##*.}
	KINGDOM=Bacteria
	# if the MAG is classified as Archaea, run Prokka with kingdom = Archaea, otherwise run with kingdom = Bacteria
	if [[ " ${ARCHAEA[@]} " =~ " ${MAG_NUM} " ]]; then
		KINGDOM=Archaea
	fi
	PROKKA_OUTDIR="$GROUP_PATH/Prokka_output/SaanichInlet_200m_${MAG_NUM}_${KINGDOM}/"
	echo "Running prokka, outputting to" $PROKKA_OUTDIR
	prokka \
	--outdir $PROKKA_OUTDIR \
	--prefix "MAG${MAG_NUM}" \
	--kingdom $KINGDOM \
	--cpus 2 \
	$MAG
	wait ;
done ;
