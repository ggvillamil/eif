


# BiocManager::install("clusterProfiler", version = "3.10")
# BiocManager::install("pathview")
# BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)


library(here)
library(magrittr)
library(tidyverse)


# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)






gtf_gr <- readRDS("../eif4f_archive/eif4f_combined2021-2022/data/gtf_gr.rds")

# table with gene 2 gname transformation
txid2gid <- gtf_gr %>%
  mcols() %>%
  as.data.frame() %>%
  distinct(transcript_id, gene_id) %>%
  select(transcript_id, gene_id) %>%
  na.omit()







subunit <- "4E"

# reading in data from deseq2
df = readRDS(here(paste0("results/post/deseq_res_deltaTE_", subunit, "_yeastnorm.rds")))






df <- df %>% as_tibble(rownames = "transcript_id") %>% left_join(txid2gid)

df <- df %>% mutate(gene_id = str_extract(gene_id, "[^.]+"))



df <- df %>% group_by(gene_id) %>% slice(which.min(padj))

# df <- df[which(df$padj < 0.05),] # hard cut-off


# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$gene_id

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)






gse <- gseGO(geneList=gene_list, 
             ont = "BP", 
             keyType = "ENSEMBL", 
             # nPerm = 10000,
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")






plotDotplot <- function(INPUT, filename){
  plotfile <- here(paste0("plots/", filename)) %T>% pdf(h=8, w=12)

  require(DOSE)
  rwplot <- dotplot(INPUT, showCategory=10, split=".sign") + facet_grid(.~.sign)

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}


gse %>% plotDotplot(., filename = "dotplot_gsea_BP_4E_lowestpadj.pdf")
