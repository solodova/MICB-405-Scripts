####################
# PATHVIEW

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pathview")
library(tidyr)
library(dplyr)
library(pathview)
library(RColorBrewer)
library(knitr)

setwd("/Users/olga/Desktop/MICB405/METAGENOME/")

ko_91 <- read.table("kaas_MAG91_ko.cleaned.txt") %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(ko = V2)

metat_rpkm_91 <- read.table("91_Bacteria_rpkm.csv", sep=",") %>%
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(rpkm = V2)

ko_201 <- read.table("kaas_MAG201_ko.cleaned.txt") %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(ko = V2)

metat_rpkm_201 <- read.table("201_Archaea_rpkm.csv", sep=",") %>%
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(rpkm = V2)

ko_169 <- read.table("kaas_MAG169_ko.cleaned.txt") %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(ko = V2)

metat_rpkm_169 <- read.table("169_Bacteria_rpkm.csv", sep=",") %>%
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(rpkm = V2)


prokka_mag_map <- read.table("Prokka_MAG_map.csv", header=F, sep=',') %>% 
  dplyr::rename(prokka_id = V1) %>% 
  dplyr::rename(mag = V2)

gtdb_dat <- rbind(arc_class, bac_class) %>% 
  separate(V2, sep=';', into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  separate(V1, sep='\\.', into=c("MAG_INFO", "MAG_ID"))

checkm_dat <- read.table("MetaBAT2_SaanichInlet_200m_min1500_checkM_stdout.tsv",
                         header=TRUE,
                         sep="\t",
                         comment.char = '') %>% 
  dplyr::rename(mag = Bin.Id) %>% 
  dplyr::select(mag, Completeness, Contamination)

# Due to a bug in the renaming script we have to rename the bins. Its a bit hacky but works using tidyverse functions
metag_rpkm <- read.table("SaanichInlet_200m_binned.rpkm.csv", header=T, sep=',') %>% 
  mutate(Sequence_name = gsub('m_', 'm.', Sequence_name)) %>% 
  mutate(Sequence_name = gsub('Inlet_', 'Inlet.', Sequence_name)) %>% 
  separate(col=Sequence_name, into=c("mag", "contig"), sep='_', extra="merge") %>% 
  group_by(Sample_name, mag) %>% 
  summarise(g_rpkm = sum(RPKM)) %>% 
  mutate(mag = gsub('Inlet.', 'Inlet_', mag))


gtdb_dat %>% 
  group_by(Phylum) %>% 
  summarise(count = n_distinct(mag)) %>% 
  kable()

gtdb_dat <- dplyr::select(gtdb_dat, mag, Kingdom, Phylum, Class, Order, Family)


rpkm_dat <- left_join(ko, metat_rpkm, by="orf") %>%
  separate(orf, into=c("prokka_id", "orf_id")) %>% # Split the Prokka ORF names into MAG identifier and ORF number for joining
  left_join(prokka_mag_map, by="prokka_id") %>% 
  left_join(gtdb_dat, by="mag") %>% 
  left_join(checkm_dat, by="mag")

# If you also wanted to add the RPKM abundance values from the metagenome:
# left_join(metag_rpkm, by="mag")

head(rpkm_dat) %>% kable()

# Nitrogen metabolism
pv.out <- pathview(gene.data = pv_mat,
                   limit = list(gene = c(0,10)),
                   low = list(gene = "#91bfdb"),
                   mid = list(gene = "#ffffbf"),
                   high = list(gene = "#fc8d59"),
                   species = "ko",
                   pathway.id="00910",
                   kegg.dir = "./")
