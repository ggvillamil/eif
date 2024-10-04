# Data preparation

library(topGO)
library(ALL)
data(ALL)
data(geneList)


affyLib <- paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)



sum(topDiffGenes(geneList))


sampleGOdata <- new("topGOdata",
  description = "Simple session", ontology = "BP",
  allGenes = geneList, geneSel = topDiffGenes,
  nodeSize = 10,
  annot = annFUN.db, affyLib = affyLib)


sampleGOdata


# Enrichment tests

resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")

resultFisher


resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")


# Analysis of results
allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
  classicKS = resultKS, elimKS = resultKS.elim, orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)


# Plot








library(topGO)
library(here)
library(tidyverse)


subunit = "4E"
dTE_down <- readRDS(here(paste0("results/post/dTE_down_", subunit, ".rds")))
dTE_down_padj <- dTE_down %>% select(gene_name, padj) %>% deframe()





# GTF GRanges copied from 2021 run (should be the same as in 2022 run)
gtf_gr <- readRDS(here("results/post/gtf_gr.rds"))


# table with gene 2 gname transformation
txid2gname <- gtf_gr %>%
  mcols() %>%
  as.data.frame() %>%
  distinct(transcript_id, gene_name) %>%
  select(transcript_id, gene_name)




# genelist <- readRDS(here("results/post/genelist.rds"))
res_dTE <- readRDS(here(paste0("results/post/deseq_res_deltaTE_", subunit, "_yeastnorm.rds")))
res_dTE <- as.data.frame(res_dTE) %>% mutate(transcript_id = rownames(.))

res_dTE <- res_dTE %>%
  left_join(txid2gname) %>%
  relocate(c(transcript_id, gene_name))

genelist <- res_dTE %>% select(gene_name, padj) %>% deframe() %>% .[complete.cases(.)]


select_dTE_down <- function(x){
  x < 0.01
}
es_dTE_down <- x %>%
    filter(padj < 0.01) %>%
    filter(log2FoldChange < 0) %>%
    arrange(log2FoldChange)


new("topGOdata",
  description = "4E_down",
  ontology = "BP",
  allGenes = genelist,
  geneSel = dTE_down_padj,
  nodeSize = 10,
  annot = annFUN.db
  )