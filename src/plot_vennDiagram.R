library(magrittr)
library(DESeq2)
library(VennDiagram)
library(tidyverse)


# Functions ---------------------------------------------------------------




# Read in differential translation results from DESeq2 --------------------


my_subunit <- "eIF3d"
my_direction <- "upregulated" # Choose "upregulated" or "downregulated"


res_files <- list.files("results/post/")
res_files <- res_files[str_detect(res_files, "deseq_res_deltaTE_yeastnorm_conditionAuxin")]
res_files <- res_files[str_detect(res_files, my_subunit)]


# Run in a loop -----


res_list <- list()
for(res_file in res_files){
  # Results from DESeq
  res_list[[res_file]] <- readRDS(paste0("results/post/", res_file)) %>%
    as_tibble(rownames = "transcript_id") %>%
    drop_na(log2FoldChange) %>%
    drop_na(padj)

  res_list[[res_file]] <- res_list[[res_file]] %>%
    mutate(diffTranslated = ifelse(log2FoldChange > 0.5 & padj < 0.05, "Up",
      ifelse(log2FoldChange < -0.5 & padj < 0.05, "Down", "No")))

  res_list[[res_file]] <- res_list[[res_file]] %>%
    select(c(transcript_id, diffTranslated))
}


# Make vectors of differential transcripts -----


up_4h <- res_list[str_detect(res_files, "4h")] %>% .[[1]] %>% filter(diffTranslated == "Up") %>% .$transcript_id
up_8h <- res_list[str_detect(res_files, "8h")] %>% .[[1]] %>% filter(diffTranslated == "Up") %>% .$transcript_id

down_4h <- res_list[str_detect(res_files, "4h")] %>% .[[1]] %>% filter(diffTranslated == "Down") %>% .$transcript_id
down_8h <- res_list[str_detect(res_files, "8h")] %>% .[[1]] %>% filter(diffTranslated == "Down") %>% .$transcript_id


# Plot Venn diagrams -----


# # Up
# x <- list(up_4h, up_8h)
# venn.diagram(x,
#   category.names = c(paste0(my_subunit, "_4h_up"), paste0(my_subunit, "_8h_up")),
#   filename = paste0("plots/venndiagram_conditionTime_", my_subunit, "_upregulated_dTE_transcripts.png"),
#   disable.logging = TRUE, imagetype = "png",
#   height = 480,
#   width = 480,
#   resolution = 300,
#   compression = "lzw",
#   lwd = 1,
#   cex = 0.3,
#   fontfamily = "sans",
#   cat.cex = 0.3,
#   cat.default.pos = "outer",
#   cat.fontfamily = "sans")


# # Down
# x <- list(down_4h, down_8h)
# venn.diagram(x,
#   category.names = c(paste0(my_subunit, "_4h_down"), paste0(my_subunit, "_8h_down")),
#   filename = paste0("plots/venndiagram_conditionTime_", my_subunit, "_downregulated_dTE_transcripts.png"),
#   disable.logging = TRUE, imagetype = "png",
#   height = 480,
#   width = 480,
#   resolution = 300,
#   compression = "lzw",
#   lwd = 1,
#   cex = 0.3,
#   fontfamily = "sans",
#   cat.cex = 0.3,
#   cat.default.pos = "outer",
#   cat.fontfamily = "sans")


# Correlate set diff log2FC of 8h vs. 4h -----


log2FC_list <- list()
for(res_file in res_files){
  # Results from DESeq
  log2FC_list[[res_file]] <- readRDS(paste0("results/post/", res_file)) %>%
    as_tibble(rownames = "transcript_id") %>%
    drop_na(log2FoldChange) %>%
    drop_na(padj)

  log2FC_list[[res_file]] <- log2FC_list[[res_file]] %>%
    select(c(transcript_id, log2FoldChange))
}


log2FC_df <- log2FC_list %>% bind_rows(.id = "time") %>%
  mutate(time = str_sub(.$time, -6, -5)) %>%
  pivot_wider(id_cols = transcript_id, names_from = time, names_prefix = "log2FoldChange_", values_from = log2FoldChange)


if(my_direction == "downregulated"){
  setdiff_tx <- union(setdiff(down_8h, down_4h), setdiff(down_4h, down_8h))
}

if(my_direction == "upregulated"){
  setdiff_tx <- union(setdiff(up_8h, up_4h), setdiff(up_4h, up_8h))
}

setdiff_log2FC <- log2FC_df %>% filter(transcript_id %in% setdiff_tx)


# Plot!
library(smplot2)

plotScatter <- function(INPUT, subtitle, filename){
  plotfile <- paste0("plots/", filename) %T>% pdf(h=8, w=8)

  # ggplot
  rwplot <- INPUT %>%
    ggplot(data = ., aes(x = log2FoldChange_4h, y = log2FoldChange_8h)) +
    geom_point() +
    sm_statCorr(corr_method = "spearman") +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    labs(title = "Correlation of log2FoldChange (8h vs. 4h)",
         subtitle = subtitle) +
    theme_bw()

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}


log2FC_df %>% plotScatter(subtitle = my_subunit, filename = paste0("scatterplot_correlation_deltaTE_log2FC_conditionAuxin_", my_subunit, "_8hvs4h.pdf"))
setdiff_log2FC %>% plotScatter(subtitle = my_subunit, filename = paste0("scatterplot_correlation_deltaTE_log2FC_conditionAuxin_", my_subunit, "_8hvs4h_setdiff_", my_direction, ".pdf"))





