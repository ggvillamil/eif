library(magrittr)
library(DESeq2)
library(VennDiagram)
library(tidyverse)


# Functions ---------------------------------------------------------------




# Read in differential translation results from DESeq2 --------------------


my_subunit <- "eIF4G3"

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


up_4h <- res_list[str_detect(res_files, "4h")] %>% .[[1]] %>% filter(diffTranslated == "Up") %>% .$transcript_id
up_8h <- res_list[str_detect(res_files, "8h")] %>% .[[1]] %>% filter(diffTranslated == "Up") %>% .$transcript_id

down_4h <- res_list[str_detect(res_files, "4h")] %>% .[[1]] %>% filter(diffTranslated == "Down") %>% .$transcript_id
down_8h <- res_list[str_detect(res_files, "8h")] %>% .[[1]] %>% filter(diffTranslated == "Down") %>% .$transcript_id


# Up
x <- list(up_4h, up_8h)
venn.diagram(x,
  category.names = c(paste0(my_subunit, "_4h_up"), paste0(my_subunit, "_8h_up")),
  filename = paste0("plots/venndiagram_conditionTime_", my_subunit, "_upregulated_dTE_transcripts.png"),
  disable.logging = TRUE, imagetype = "png",
  height = 480,
  width = 480,
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  cex = 0.3,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.fontfamily = "sans")


# Down
x <- list(down_4h, down_8h)
venn.diagram(x,
  category.names = c(paste0(my_subunit, "_4h_down"), paste0(my_subunit, "_8h_down")),
  filename = paste0("plots/venndiagram_conditionTime_", my_subunit, "_downregulated_dTE_transcripts.png"),
  disable.logging = TRUE, imagetype = "png",
  height = 480,
  width = 480,
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  cex = 0.3,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.fontfamily = "sans")
