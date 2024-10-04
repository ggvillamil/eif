

library(here)
library(magrittr)
library(LSD)
library(tidyverse)



featurestbl <- readRDS(here("results/post/featurestbl.rds"))



featurestbl <- featurestbl %>%
  mutate(transcript_id = transcript_id_version) %>%
  dplyr::select(-c(transcript_id_version))


dTE_down <- readRDS(here(paste0("results/post/top500_dTE_down_4G2.rds"))) %>%
  select(c(transcript_id, gene_name, log2FoldChange, padj))

dTE_up <- readRDS(here(paste0("results/post/top500_dTE_up_4G2.rds"))) %>%
  select(c(transcript_id, gene_name, log2FoldChange, padj))


dTE <- rbind(dTE_up[1:100,], dTE_down[1:100,])


dTE_features <- dTE %>%
  left_join(featurestbl)

dTE_down_features <- dTE_down %>%
  left_join(featurestbl)



plotHeatscatter <- function(INPUT, filename){
  plotfile <- here(paste0("plots/", filename)) %T>% pdf(h=6, w=6)

  rwplot <- heatscatter(x = INPUT$gc_fputrs, y = INPUT$log2FoldChange, cor = TRUE)

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}



dTE_down_features %>% plotHeatscatter(., filename = "test2.pdf")