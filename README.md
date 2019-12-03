# MICB 405 Scripts
## Project 2
This folder contains scripts used for metagenomic analysis. 
- `plot_cont_vs_comp.R` can be used to build contamination vs completion plots from MAG checkM output, rpkm values, and gtdb classifications for MAGs.
- `gunzip.sh` is used to unzip metatranscriptomic .fastq files
- `prokka_200m.sh` is used to obtain prokka-predicted ORFs given .fa file corresponding to particular MAG
- `bwa.sh` is used to generate index files for given MAG (using it's .ffn prokka output file) and then run bwa mem to align metatranscriptomic reads to the reference
- `rpkm.sh` is used to generate transcriptomic rpkm values for a given MAG, for each of its prokka-predicted ORFs
- `pathview.r` is used to call pathview with data containing rpkm values from prokka ORFs mapped to KO numbers from KAAS, allowing visualization of MAG transcriptomic activity for specific KEGG pathways
- `plot_geochemical_gradients.r` is used to plot geochemical gradients
