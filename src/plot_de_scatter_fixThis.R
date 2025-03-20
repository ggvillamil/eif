library(magrittr)
library(DESeq2)
library(tidyverse)


# Read in differential translation results from DESeq2 --------------------


res_files <- list.files("results/post/")
res_files <- res_files[str_detect(res_files, "deseq_res_")]




res_ribo <- readRDS("results/post/deseq_res_diffribo_yeastnorm_conditionAuxin_human_morf_eIF4G1_minusHarr_4h.rds") %>%
  as_tibble(rownames = "transcript_id") %>%
  drop_na(log2FoldChange) %>%
  drop_na(padj) %>%
  select(c(transcript_id, log2FoldChange, padj)) %>%
  rename(l2FC_ribo = "log2FoldChange") %>%
  rename(padj_ribo = "padj")


res_total <- readRDS("results/post/deseq_res_difftotal_autonorm_conditionAuxin_human_morf_eIF4G1_minusHarr_4h.rds") %>%
  as_tibble(rownames = "transcript_id") %>%
  drop_na(log2FoldChange) %>%
  drop_na(padj) %>%
  select(c(transcript_id, log2FoldChange, padj)) %>%
  rename(l2FC_total = "log2FoldChange") %>%
  rename(padj_total = "padj")


res_combined <- left_join(res_ribo, res_total) %>%
  mutate(label = ifelse(l2FC_ribo > 0.5 & padj_ribo < 0.05, "Up",
  	ifelse(l2FC_ribo < -0.5 & padj_ribo < 0.05, "Down", "No Change"))) %>%
  mutate(label = as.factor(label))




# Function
plotScatter <- function(INPUT, subtitle, filename, xlim, ylim){
  plotfile <- paste0("plots/", filename) %T>% pdf(h=6, w=6)

  # ggplot
  rwplot <- INPUT %>%
    ggplot(data = ., aes(x = l2FC_total, y = l2FC_ribo, color = label)) +
    geom_point() +
    scale_color_manual(values = c("red", "gray", "blue")) +
    scale_alpha_manual(values = c(1, 0.5, 1)) +
    labs(title = "Change in ribosome occupancy vs. change in gene expression",
         subtitle = subtitle) +
    ylab("log2FoldChange Ribosome Occupancy") +
    xlab("log2FoldChange Gene Expression") +
    xlim(xlim) +
    ylim(ylim) +
    theme_bw()

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}


res_combined %>% plotScatter(., subtitle = "eIF4G1", filename = "scatter_log2FCriboVStotal_conditionAuxin_human_morf_eIF4G1_minusHarr_4h.pdf", xlim = c(-10,10), ylim = c(-10,10))




