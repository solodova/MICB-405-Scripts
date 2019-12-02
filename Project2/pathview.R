if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pathview")
library(tidyr)
library(dplyr)
library(pathview)
library(RColorBrewer)
library(knitr)

# set working directory to location of data
setwd("/Users/olga/Desktop/MICB405/METAGENOME/")

# read in KAAS output (with prokka ORFs that don't map to KO numbers removed),
# and metatranscriptomic rpkms found using prokka-predicted ORFs and .fastq read files
# for the MAGs that we have decided to pursue (the 3 HQ MAGs)

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

ko <- rbind(ko_91, ko_169, ko_201)
metat_rpkm <- rbind(metat_rpkm_91, metat_rpkm_169, metat_rpkm_201)

# read in table mapping prokka IDs to MAG IDs
prokka_mag_map <- read.table("Prokka_MAG_map.csv", header=F, sep=',') %>% 
  dplyr::rename(prokka_id = V1) %>% 
  dplyr::rename(MAG_ID = V2) 
prokka_mag_map$MAG_ID <- as.character(prokka_mag_map$MAG_ID)

# join rpkm data with KO numbers using prokka ORF ID
rpkm_dat <- left_join(ko, metat_rpkm, by="orf") %>%
  separate(orf, into=c("prokka_id", "orf_id")) %>% 
  left_join(prokka_mag_map, by="prokka_id") 

# to look at output for specific MAG, uncomment filter operator
# pv_mat contains t_rpkm values for each KO number associated with a MAG
pv_mat <- rpkm_dat %>% 
  group_by(MAG_ID, ko) %>% 
  #filter(MAG_ID == "91") %>%
  summarise(t_rpkm = sum(rpkm)) %>% 
  spread(key = MAG_ID, value = t_rpkm) %>%
  column_to_rownames(var="ko")

# Methane metabolism: 00680
# Nitrogen metabolism: 00910
# Sulphur metabolism: 00920
# call pathview using some pathway ID and rpkm values associated with KO numbers
pv.out <- pathview(gene.data = pv_mat,
                   limit = list(gene = c(0,10)),
                   low = list(gene = "#91bfdb"),
                   mid = list(gene = "#ffffbf"),
                   high = list(gene = "#fc8d59"),
                   species = "ko",
                   pathway.id="00920",
                   kegg.dir = "./pathview_output/")
