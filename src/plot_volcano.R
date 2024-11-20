library(magrittr)
library(DESeq2)
library(tidyverse)


# Functions ---------------------------------------------------------------


plotVolcano <- function(INPUT, subtitle, filename, xlim, ylim){
  plotfile <- paste0("plots/", filename) %T>% pdf(h=8, w=14)

  # ggplot
  rwplot <- INPUT %>%
    ggplot(data = ., aes(x = log2FoldChange, y = -log10(padj), color = diffTranslated)) +
    geom_point() +
    scale_color_manual(values = c("red", "black", "blue")) +
    labs(title = "Differential Translational Efficiency",
         subtitle = subtitle) +
    ylab("-Log10 (adjusted p-value)") +
    xlab("TE log2FoldChange") +
    xlim(xlim) +
    ylim(ylim) +
    theme_bw()

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}


# Read in differential translation results from DESeq2 --------------------


res_files <- list.files("results/post/")
res_files <- res_files[str_detect(res_files, "deseq_res_")]


# Run in a loop -----


for(res_file in res_files){
  subtitle <- str_sub(res_file, 11, -5)
  filename <- paste0("volcanoplot_", str_sub(res_file, 1, -5), ".pdf")

  # Results from DESeq
  res <- readRDS(paste0("results/post/", res_file)) %>%
    as_tibble(rownames = "transcript_id") %>%
    drop_na(log2FoldChange) %>%
    drop_na(padj)

  res <- res %>%
    mutate(diffTranslated = ifelse(log2FoldChange > 0.5 & padj < 0.05, "Up",
      ifelse(log2FoldChange < -0.5 & padj < 0.05, "Down", "No")))

  res %>% plotVolcano(., subtitle = subtitle, filename = filename, xlim = c(-15, 15), ylim = c(-1, 30))
}














