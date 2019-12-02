library(tidyverse)

# set working directory to folder containing data
setwd("/Users/olga/Desktop/MICB405/METAGENOME/")

# read in classifications from gtdbtk for medium and high quality MAGs 
ar122_class <- read.csv("gtdbtk.ar122.classification_pplacer.tsv", 
                        sep="\t", 
                        header=FALSE)
colnames(ar122_class) <- c("Bin.Id", "classification")
bac120_class <- read.csv("gtdbtk.bac120.classification_pplacer.tsv", 
                         sep="\t", 
                         header=FALSE)
colnames(bac120_class) <- c("Bin.Id", "classification")
# read in checkM output with quality check (contamination/completion) for all MAGs
checkM_output <- read.csv("MetaBAT2_SaanichInlet_200m_min1500_checkM_stdout.tsv", 
                         sep="\t")

MedQHigh_MAGs <- c(as.character(ar122_class$Bin.Id), 
                   as.character(bac120_class$Bin.Id))

# keep checkM output only for medium/high quality MAGS
# then join checkM output with classifications for each MAG
checkM_gtdbtk_joined <- checkM_output %>% 
                      filter(Bin.Id %in% MedQHigh_MAGs) %>%
                      left_join(ar122_class) %>% 
                      left_join(bac120_class, by="Bin.Id")
checkM_gtdbtk_joined$classification.x <- as.character(checkM_gtdbtk_joined$classification.x)
checkM_gtdbtk_joined$classification.y <- as.character(checkM_gtdbtk_joined$classification.y)

# set classification of a MAG (classification.x) to be either bacteria/archaea classification from gtdbtk
# (whichever one was not null for the given MAG)
checkM_gtdbtk_joined$classification.x[is.na(checkM_gtdbtk_joined$classification.x)] <- 
  checkM_gtdbtk_joined$classification.y[is.na(checkM_gtdbtk_joined$classification.x)]
checkM_gtdbtk_joined$classification.y <- NULL

# separate columns for classification and MAG identifiers
checkM_gtdbtk_joined <- checkM_gtdbtk_joined %>% 
  separate(col="classification.x", into=c("d", "p", "c", "o", "f", "g", "s"), sep=";") %>%
  separate(col="Bin.Id", into=c("Depth/Area", "MAG_ID"), sep="\\.") 

# clean up the classifications in the table
for (col in c("d", "p", "c", "o", "f", "g", "s")) {
  checkM_gtdbtk_joined[, col] <- str_remove(checkM_gtdbtk_joined[, col], ".__")
}

# set best_class (best classification) to be the most specific taxonomy available for each MAG
checkM_gtdbtk_joined$best_class <- c("Anaerolineales",       
                                     "Gracilibacteria",
                                     "Komeilibacterales",
                                     "Pacearchaeales",
                                     "Elusimicrobiales",
                                     "Pseudomonadales",
                                     "Alphaproteobacteria",
                                     "Bacteroidota",
                                     "Prolixibacteraceae",
                                     "Komeilibacterales",
                                     "Thioglobaceae",
                                     "Nitrosopelagicus",
                                     "Magasanikibacterales",
                                     "Gammaproteobacteria",
                                     "Woesearchaeia",
                                     "Woesearchaeia",
                                     "Anaerolineales",
                                     "Flavobacteriales",
                                     "Pseudohongiellaceae",
                                     "Woesearchaeia",
                                     "Methyloprofundus",
                                     "Nitrosomonadaceae",
                                     "Marinisomatales",
                                     "Moranbacterales",
                                     "Pacearchaeales",
                                     "Bathyarchaeia",
                                     "Gracilibacteria",
                                     "Iainarchaeales",
                                     "Bacteria",
                                     "Alphaproteobacteria",
                                     "Pacebacteraceae",
                                     "Woesearchaeia",
                                     "Komeilibacterales",
                                     "Babeliales",
                                     "Pseudomonadales",
                                     "Cryomorphaceae",
                                     "Nitrosopumilus",
                                     "Woesearchaeia",
                                     "Nanoarchaeota",
                                     "Pacearchaeales",
                                     "Berkiellales",
                                     "Woesearchaeia",
                                     "Woesearchaeia",
                                     "Amylibacter",
                                     "Woesearchaeia",
                                     "Peregrinibacterales",
                                     "Gammaproteobacteria",
                                     "Bacteroidales",
                                     "Sneathiellales",
                                     "Gemmatimonadetes")
checkM_gtdbtk_sorted_taxonomy <- 
  checkM_gtdbtk_joined[order(checkM_gtdbtk_joined[,"d"], 
                             checkM_gtdbtk_joined[,"p"], 
                             checkM_gtdbtk_joined[,"c"],
                             checkM_gtdbtk_joined[,"o"],
                             checkM_gtdbtk_joined[,"f"],
                             checkM_gtdbtk_joined[,"g"],
                             checkM_gtdbtk_joined[,"s"]),]
  

# get the genomic rpkms for MAGs/contigs, sum rpkm values over given MAG
genomic_rpkms_SI072 <- read.csv("SaanichInlet_200m_binned.rpkm.csv")
genomic_rpkms_SI072 <- genomic_rpkms_SI072 %>% 
  separate(col="Sequence_name", into=c("Area", "Depth", "MAG_ID", "Contig_ID"), sep="_") 
genomic_rpkms_summary <- aggregate(RPKM ~ MAG_ID, data=genomic_rpkms_SI072, FUN=sum)

# get transcriptomic rpkms for MAGs/contigs, sum rpkm values over given MAG
transcriptomic_rpkms_SI072 <- read.csv("SI072_200m_metaT_rpkm.csv")
transcriptomic_rpkms_SI072 <- transcriptomic_rpkms_SI072 %>% 
  separate(col="Sequence", into=c("Area", "Depth", "MAG_ID", "Contig_ID"), sep="_") 
transcriptomic_rpkms_summary <- aggregate(RPKM ~ MAG_ID, data=transcriptomic_rpkms_SI072, FUN=sum)

# add rpkm info to table containing checkM and gtdbtk outputs
checkM_gtdbtk_joined<- checkM_gtdbtk_sorted_taxonomy %>% 
  left_join(genomic_rpkms_summary) %>% 
  left_join(transcriptomic_rpkms_summary, by="MAG_ID")
colnames(checkM_gtdbtk_joined) <- 
  c(colnames(checkM_gtdbtk_joined)[1:(length(colnames(checkM_gtdbtk_joined)) - 2)], "RPKM_genomic", "RPKM_transcriptomic")
write.csv(checkM_gtdbtk_joined, "checkM_gtdbtk_rpkm_200mMAGs.csv")

# function which will create contamination vs completion plot, each MAG as a point with rpkm 
# (genomic or transcriptomic, whichever specified) encoded as size
plot_MAGs <- function(checkM_gtdbtk, classifications, taxonomy, rpkm, rpkm_type) {
  plot <- checkM_gtdbtk %>%
    ggplot() +
    geom_point(aes(x=Completeness, 
                   y=Contamination,
                   color=classifications, 
                   shape=d, size=rpkm)) +
    scale_size_continuous(range = c(1,4)) +
    labs(title ="Completeness vs Contamination for Medium/High Quality MAGs from sample SI072 at 200m", 
         x = "Completeness (%)", y = "Contamination (%)",
         color=paste("Taxonomy (",taxonomy ,")",sep=""), 
         shape="Domain",
         size=paste(rpkm_type, " RPKM")) +
    annotate("rect", xmin=90, xmax=100, ymin=0, ymax= 5, 
                 fill=NA, colour="red") +
    #ylim(0, 100) +
    theme_bw() +
    theme(legend.position = "bottom") +
    theme(legend.direction = "horizontal") +
    theme(legend.title.align = 0.5) +
    theme(legend.box = "vertical") +
    guides(shape=guide_legend(order=1),
           size=guide_legend(order=2),
           color=guide_legend(nrow=5, order=3, title.position = "top")) +
    theme(legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")) +
    theme(legend.key.size = unit(0.5, "cm"))
    
  return(plot)
}

# create contamination vs completion plots using both genomic and transcriptomic rpkms,
# also using either phyum or best classification for MAG as it's taxonomy
phylum_genomic <- plot_MAGs(checkM_gtdbtk=checkM_gtdbtk_joined,
                            classifications=checkM_gtdbtk_joined$p, 
                            taxonomy="phylum", 
                            rpkm=checkM_gtdbtk_joined$RPKM_genomic, 
                            rpkm_type="Genomic")
phylum_transcriptomic <-plot_MAGs(checkM_gtdbtk=checkM_gtdbtk_joined, 
                                  classifications=checkM_gtdbtk_joined$p, 
                                  taxonomy="phylum", 
                                  rpkm=checkM_gtdbtk_joined$RPKM_transcriptomic, 
                                  rpkm_type="Transcriptomic")

bestclass_genomic <- plot_MAGs(checkM_gtdbtk=checkM_gtdbtk_joined,
                            classifications=checkM_gtdbtk_joined$best_class, 
                            taxonomy="most specific available", 
                            rpkm=checkM_gtdbtk_joined$RPKM_genomic, 
                            rpkm_type="Genomic")
bestclass_transcriptomic <-plot_MAGs(checkM_gtdbtk=checkM_gtdbtk_joined, 
                                  classifications=checkM_gtdbtk_joined$best_class, 
                                  taxonomy="most specific available", 
                                  rpkm=checkM_gtdbtk_joined$RPKM_transcriptomic, 
                                  rpkm_type="Transcriptomic")

# save just the high quality MAGs
HQ_MAGs <- checkM_gtdbtk_joined %>% filter((Contamination < 5) & (Completeness > 90))
write.csv(HQ_MAGs, "checkM_gtdbtk_rpkm_200mHQMAGs.csv")

# for both just high quality, and for high + medium quality, make bar plots showing transcriptomic and
# genomic rpkms for each MAG
HQ_MAG_rpkm <- HQ_MAGs[order(as.integer(HQ_MAGs$MAG_ID)),] %>%
      gather(key="RPKM", value="Value", RPKM_genomic, RPKM_transcriptomic) %>% 
      ggplot(aes(x=MAG_ID,y=Value,fill=factor(RPKM)))+
      geom_bar(stat="identity",position="dodge")+
      labs(title ="RPKM for High Quality MAGs from sample SI072 at 200m", 
                                                     x = "MAG ID", y = "RPKM value",
                                                     fill="RPKM type")
checkM_gtdbtk_joined <- checkM_gtdbtk_joined[order(as.integer(checkM_gtdbtk_joined$MAG_ID)),] %>%
  gather(key="RPKM", value="Value", RPKM_genomic, RPKM_transcriptomic) %>% 
  ggplot(aes(x=MAG_ID,y=Value,fill=factor(RPKM)))+
  geom_bar(stat="identity",position="dodge")+
  labs(title ="RPKM for High Quality MAGs from sample SI072 at 200m", 
       x = "MAG ID", y = "RPKM value",
       fill="RPKM type")


