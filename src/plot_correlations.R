library(magrittr)
library(DESeq2)
library(tidyverse)
library(LSD)


# Read in differential translation results from DESeq2 --------------------


res_files <- list.files("results/post/")
res_files <- res_files[str_detect(res_files, "deseq_res_")]




res_eif3d <- readRDS("results/post/deseq_res_deltaTE_yeastnorm_conditionAuxin_human_morf_eIF3d_minusHarr_4h.rds") %>%
  as_tibble(rownames = "transcript_id") %>%
  drop_na(log2FoldChange) %>%
  drop_na(padj) %>%
  select(c(transcript_id, log2FoldChange, padj)) %>%
  rename(l2FC_eif3d = "log2FoldChange") %>%
  rename(padj_eif3d = "padj")


res_eif4e <- readRDS("results/post/deseq_res_deltaTE_yeastnorm_conditionAuxin_human_morf_eIF4E_minusHarr_4h.rds") %>%
  as_tibble(rownames = "transcript_id") %>%
  drop_na(log2FoldChange) %>%
  drop_na(padj) %>%
  select(c(transcript_id, log2FoldChange, padj)) %>%
  rename(l2FC_eif4e = "log2FoldChange") %>%
  rename(padj_eif4e = "padj")


res_eif4g1 <- readRDS("results/post/deseq_res_deltaTE_yeastnorm_conditionAuxin_human_morf_eIF4G1_minusHarr_4h.rds") %>%
  as_tibble(rownames = "transcript_id") %>%
  drop_na(log2FoldChange) %>%
  drop_na(padj) %>%
  select(c(transcript_id, log2FoldChange, padj)) %>%
  rename(l2FC_eif4g1 = "log2FoldChange") %>%
  rename(padj_eif4g1 = "padj")


res_3d4e <- left_join(res_eif3d, res_eif4e) %>% filter(padj_eif3d < 0.01 | padj_eif4e < 0.01)
res_3d4g1 <- left_join(res_eif3d, res_eif4g1) %>% filter(padj_eif3d < 0.01 | padj_eif4g1 < 0.01)
res_4e4g1 <- left_join(res_eif4e, res_eif4g1) %>% filter(padj_eif4e < 0.01 | padj_eif4g1 < 0.01)



# Function
plotScatter <- function(x, y, filename){
  plotfile <- paste0("plots/", filename) %T>% pdf(h=6, w=6)

  rwplot <- heatscatter(x, y, cor = TRUE, method = "pearson",
    xlab = "log2FoldChange TE eIF4E",
    ylab = "log2FoldChange TE eIF4G1",
    # xlim = c(-10,10),
    # ylim = c(-10,10),
    main = "Correlation of TE fold change (log2)")

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}


plotScatter(x = res_4e4g1$l2FC_eif4e, y = res_4e4g1$l2FC_eif4g1, filename = "scatter_log2FC_correlation_eIF4E_vs_eIF4G1.pdf")




