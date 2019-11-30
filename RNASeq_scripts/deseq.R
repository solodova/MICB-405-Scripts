if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("biomaRt")
BiocManager::install("ReactomePA")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(magrittr)
library(org.Mm.eg.db)
library(ReactomePA)

mart <- biomaRt::useDataset("mmusculus_gene_ensembl", biomaRt::useMart("ensembl"))

##### function to map gene ids to descriptions of the genes
mapStuff <- function(genes, ID) {
  if (ID == "ensembl"){
    geneList <- biomaRt::getBM(filters= "ensembl_gene_id", 
                               attributes= c("ensembl_gene_id", "description"),
                               values= as.character(genes),
                               mart=mart)
    # temp contains genes which are in genes but not in geneList
    temp <- data.frame(ensembl_gene_id= setdiff(as.character(genes), geneList$ensemble_gene_id))
    geneList %<>% full_join(temp)
    rm(temp)
    geneList <- geneList[match(as.character(genes), geneList$ensembl_gene_id),]
    return(geneList)
  }
}

##### set up directory
setwd("/Users/olga/Desktop/MICB405/")

##### set meta data for counts
sampleNames <- c("PRIMARY1", "PRIMARY2", "PRIMARY3", "NAIVE1", "NAIVE2", "NAIVE3")
conditions <- c("primary", "primary", "primary", "naive", "naive", "naive")
fileNames <- paste0(paste0("htseq_counts_", sampleNames), ".txt")
sample_metadata <- data.frame(sampleName=sampleNames,
                              fileName=fileNames,
                              condition=conditions)

##### make DESeqDataSet and do DESeq
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sample_metadata,
                                       design = ~ condition)
ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "naive")
dds <- DESeq(ddsHTSeq)

##### DEseq results
res <- results(dds, name = "condition_primary_vs_naive")

write.csv(as.data.frame(subset(res, padj < 0.05)), 
          file = "primary_vs_naive_05_vanillaDEoutput.txt")

##### plot counts for particular gene across the dataset
plotCounts(dds, gene="ENSMUSG00000034390", intgroup="condition", main="c-Maf")
plotCounts(dds, gene="ENSMUSG00000021356", intgroup="condition", main="IRF4")

##### create results dataframe for genes which whose padj < 0.05, add additional columns
PADJ <- 0.05
significant_res <- 
  res %>% as.data.frame() %>% 
  rownames_to_column(var= "ensembl_gene_id") %>% 
  mutate(Abs_LFC = abs(log2FoldChange)) %>% 
  mutate(fc= sign(log2FoldChange)* 2^(Abs_LFC)) %>% 
  subset((padj < PADJ))

rownames(significant_res) <- NULL
save_res <- significant_res
save_res <- column_to_rownames(save_res, var="ensembl_gene_id")
save_res <- subset(save_res, padj < 0.001)
save_res_n <- subset(save_res, fc < -2)
save_res_p <- subset(save_res, fc > 2)

save_res_n <- save_res_n[, !names(save_res_n) %in% c("Abs_LFC", "fc")]
save_res_p <- save_res_p[, !names(save_res_p) %in% c("Abs_LFC", "fc")]

write.csv(as.data.frame(save_res_p),
          file = "primary_vs_naive_001_fc2.csv")
write.csv(as.data.frame(save_res_n),
          file = "primary_vs_naive_001_fc-2.csv")
rm(save_res)
rm(save_res_n)
rm(save_res_p)
##### add gene descriptions to results dataframe
significant_res <- 
  significant_res$ensembl_gene_id %>% 
  mapStuff(ID = "ensembl") %>% 
  right_join(significant_res)

##### all ensembl to entrex conversion
##### 16% of genes fail to map here
universe_dups <- 
  significant_res$ensembl_gene_id %>% 
  clusterProfiler::bitr(fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
colnames(universe_dups) <- c("ensembl_gene_id", "entrez_gene_id")
significant_res_dup <- universe_dups %>% right_join(significant_res)

up <- subset(significant_res_dup, fc >=1.5)
down <- subset(significant_res_dup, fc <= -1.5)

x <- enrichPathway(gene = up$entrez_gene_id[!is.na(up$entrez_gene_id)], 
                   universe = as.character(universe_dups$entrez_gene_id), 
                   organism="mouse", pvalueCutoff=0.05, readable=T)

dotplot(x, showCategory=15)


##### remove duplicated entries (ensembl ids with more than one entrez id)
dup_ids <- subset(data.frame(table(universe_dups$ensembl_gene_id)), Freq > 1)
universe_nodups <- subset(universe_dups, !(ensembl_gene_id %in% dup_ids$Var1))
significant_res_nodup <- universe_nodups %>% right_join(significant_res)

##### save csv file for all genes, only upregulated, only downregulated
##### order results by absolute fold change
significant_res_save <- significant_res_dup[order(-(abs(significant_res_dup$fc))),]
write.csv(as.data.frame(significant_res_save), 
          file = "primary_vs_naive_05.csv")
write.csv(as.data.frame(subset(significant_res_save, fc > 0)), 
          file = "primary_vs_naive_05_upr.csv")
write.csv(as.data.frame(subset(significant_res_save, fc < 0)), 
          file = "primary_vs_naive_05_downr.csv")


##### PATHWAY ENRICHMENT

paths_dup <- get_dotplot_enriched_paths(significant_res_dup, up=TRUE, down=TRUE, universe=universe_dups)
paths_nodup <- get_dotplot_enriched_paths(significant_res_nodup, up=TRUE, down=TRUE, universe=universe_nodups)

paths_dup_upr <- get_dotplot_enriched_paths(significant_res_dup, up=TRUE, down=FALSE, universe=universe_dups)
write.csv(as.data.frame(paths_dup_upr$signif_paths), 
          file = "primary_vs_naive_05_upr_reactomedup.csv")
paths_nodup_upr <- get_dotplot_enriched_paths(significant_res_nodup, up=TRUE, down=FALSE, universe=universe_nodups)
write.csv(as.data.frame(paths_dup_upr$signif_paths), 
          file = "primary_vs_naive_05_upr_reactomenodup.csv")
paths_dup_downr <- get_dotplot_enriched_paths(significant_res_dup, up=FALSE, down=TRUE, universe=universe_dups)
write.csv(as.data.frame(paths_dup_upr$signif_paths), 
          file = "primary_vs_naive_05_downr_reactomedup.csv")
paths_nodup_downr <- get_dotplot_enriched_paths(significant_res_nodup, up=FALSE, down=TRUE, universe=universe_nodups)
write.csv(as.data.frame(paths_dup_upr$signif_paths), 
          file = "primary_vs_naive_05_downr_reactomenodup.csv")

get_dotplot_enriched_paths <- function(DE_signif_res, FC=1.5, PADJ = 0.05, pathPADJ = 0.05, up, down, universe) {
  DE_signif <- DE_signif_res
  if (up & down) { DE_signif <- DE_signif %>% subset((fc >= FC) | (fc <= -(FC)))}
  if (up & !down) { DE_signif <- DE_signif %>% subset(fc >= FC)}
  if (down & !up) { DE_signif <- DE_signif %>% subset(fc <= (-FC))}
  
  all_paths <- enrichPathway(gene = DE_signif$entrez_gene_id[!is.na(DE_signif$entrez_gene_id)], 
                       universe = as.character(universe$entrez_gene_id), 
                       organism="mouse")@result
  signif_paths <- all_paths[all_paths$p.adjust <= pathPADJ,]
  
  dot_plot <- dplyr::select(signif_paths, one_of("Description", "p.adjust"))
  dot_plot$Description <- factor(dot_plot$Description, levels = (dot_plot$Description %>% table() %>% sort(decreasing = FALSE) %>% names())) 
  dot_plot <- dot_plot %>% ggplot(aes(x=p.adjust, y = Description)) + geom_point() 
  
  return(list(dot_plot=dot_plot, signif_paths=signif_paths))
}





##### data visualization
rld <- rlog(dds)
plotPCA(rld, intgroup = "condition")

sample_dists <- dist(t(assay(rld)))
sample_dist_matrix <- as.matrix(sample_dists)
colnames(sample_dist_matrix) <- NULL
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists, 
         col = colours)
