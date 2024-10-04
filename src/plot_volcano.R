library(here)
library(abind)
library(magrittr)
library(ggpubr)
library(ggrepel)
library(rstatix)
library(DESeq2)
library(tidyverse)




#####


# GTF GRanges copied from 2021 run (should be the same as in 2022 run)
gtf_gr <- readRDS(here("results/post/gtf_gr.rds"))


# # Table with gene ID to gene name transformation
# gid2gname <- gtf_gr %>%
#   mcols %>%
#   as.data.frame %>%
#   distinct(gene_id, gene_name) %>%
#   select(gene_id, gene_name)

# table with gene 2 gname transformation
txid2gname <- gtf_gr %>%
  mcols() %>%
  as.data.frame() %>%
  distinct(transcript_id, gene_name) %>%
  select(transcript_id, gene_name)


saveRDS(txid2gname, "results/post/txid2gname.rds")


# Read in differential translation results from DESeq2 --------------------


txid2gname <- readRDS("results/post/txid2gname.rds")

subunit <- "4G3"
# subunit <- "oldmonosome4G1"

res_file <- paste0("results/post/deseq_res_deltaTE_", subunit, "_yeastnorm.rds")

# Results from DESeq
res <- readRDS(here(res_file)) %>%
  as_tibble(rownames = "transcript_id") %>%
  drop_na(log2FoldChange) %>%
  drop_na(padj)






res <- res %>%
  left_join(txid2gname)


#####


res <- res %>%
  mutate(diffTranslated = ifelse(log2FoldChange > 0.5 & padj < 0.05, "Up",
    ifelse(log2FoldChange < -0.5 & padj < 0.05, "Down", "No"))) %>%
  mutate(diffLabel = ifelse(diffTranslated == "No", NA, gene_name))


#####


my_genes <- readRDS(here("results/post/genelist_initiationfactors.rds"))
my_genes <- c(my_genes,"ATF4") # Integrated stress response gene ATF4 (transcription factor)
my_genes <- ""

res <- res %>%
  mutate(my_genes = ifelse(gene_name %in% my_genes, gene_name, NA))


# res <- res %>% mutate(my_genes = ifelse(diffTranslated == "No" & my_genes != NA, NA, gene_name))


#####



most_abundant_txs <- read_tsv(here(paste0("results/post/most_abundant_txs_", subunit, ".tsv")))

my_transcripts <- most_abundant_txs$transcript_id

res <- res %>% filter(transcript_id %in% my_transcripts)


#####




plotVolcano <- function(INPUT, subtitle, filename, xlim, ylim){
  plotfile <- here(paste0("plots/", filename)) %T>% pdf(h=8, w=14)

  # ggplot
  rwplot <- INPUT %>%
    ggplot(data = ., aes(x = log2FoldChange, y = -log10(padj), color = diffTranslated, label = my_genes)) +
    geom_point() +
    geom_label() +
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


res %>%
  plotVolcano(., subtitle = subunit, filename = paste0("volcanoplot_log2FCTE_mostAbundantTx_", subunit, "_initiationfactors.pdf"), xlim = c(-11, 11), ylim = c(-1, 25))



res <- res %>% mutate(my_genes = ifelse(diffTranslated != "No", gene_name, NA))


res %>%
  plotVolcano(., subtitle = subunit, filename = paste0("volcanoplot_log2FCTE_mostAbundantTx_", subunit, "_upAndDown.pdf"), xlim = c(-11, 11), ylim = c(-1, 25))





# #####


# res <- res %>%
#   mutate(membrane_specific = ifelse(gene_name %in% my_genes, "membrane_specific", NA))




# plotVolcano <- function(INPUT, subtitle, filename, xlim, ylim){
#   plotfile <- here(paste0("plots/", filename)) %T>% pdf(h=8, w=14)

#   # ggplot
#   rwplot <- INPUT %>%
#     ggplot(data = ., aes(x = log2FoldChange, y = -log10(padj), color = membrane_specific)) +
#     geom_point(alpha = 0.5) +
#     scale_color_manual(values = c("red", "black")) +
#     labs(title = "Differential Translational Efficiency",
#          subtitle = subtitle) +
#     ylab("-Log10 (adjusted p-value)") +
#     xlab("TE log2FoldChange") +
#     xlim(xlim) +
#     ylim(ylim) +
#     theme_bw()

#   print(rwplot)
#   dev.off()
#   normalizePath(plotfile) %>% message
# }


# res %>%
#   plotVolcano(., subtitle = "oldmonosome4G1", filename = "volcanoplot_log2FCTE_allGenes_oldmonosome4G1_membranespecific.pdf", xlim = c(-20, 20), ylim = c(-1, 26))







