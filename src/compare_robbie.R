library(here)
library(tidyverse)
library(GenomicRanges)
library(LSD)
library(magrittr)



res_prot <- read.csv(here("from_Robbie/g2_nsp.csv"), header = TRUE)

res_prot <- res_prot %>%
  mutate(Protein.Group = str_extract(Protein.Group, "\\w+$")) %>%
  dplyr::rename(gene_name = Protein.Group)

log2fc_prot <- res_prot %>%
  # filter(pval < 0.05) %>%
  select(gene_name, log2fc)






# GTF GRanges copied from 2021 run (should be the same as in 2022 run)
gtf_gr <- readRDS(here("results/post/gtf_gr.rds"))


# table with gene 2 gname transformation
txid2gname <- gtf_gr %>%
  mcols() %>%
  as.data.frame() %>%
  distinct(transcript_id, gene_name) %>%
  select(transcript_id, gene_name)







res_dTE <- readRDS(here("results/post/deseq_res_deltaTE_4G2_yeastnorm.rds"))
res_dTE <- as.data.frame(res_dTE) %>% mutate(transcript_id = rownames(.))


res_dTE <- res_dTE %>%
  left_join(txid2gname) %>%
  relocate(c(transcript_id, gene_name))


log2fc_dTE <- res_dTE %>%
  select(gene_name, log2FoldChange, padj) %>%
  dplyr::rename(log2fc = log2FoldChange)


# log2fc_dTE <- log2fc_dTE %>%
#   filter(padj < 0.05)


common <- intersect(log2fc_prot$gene_name, log2fc_dTE$gene_name)

log2fc_dTE <- log2fc_dTE %>%
  filter(gene_name %in% common)


# Select lowest padj transcript
log2fc_dTE <- log2fc_dTE %>% group_by(gene_name) %>% slice_min(order_by = padj, n = 1) %>% select(-padj)



log2fc_prot <- log2fc_prot %>% filter(gene_name %in% common)



log2fc_dTE <- log2fc_dTE %>% dplyr::rename(log2fc_dTE = log2fc)
log2fc_prot <- log2fc_prot %>% dplyr::rename(log2fc_prot = log2fc)


INPUT <- log2fc_dTE %>% left_join(log2fc_prot)


plotHeatscatter <- function(INPUT, filename){
  plotfile <- here(paste0("plots/", filename)) %T>% pdf(h=6, w=6)

  rwplot <- heatscatter(x = INPUT$log2fc_prot, y = INPUT$log2fc_dTE, cor = TRUE)

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}



INPUT %>% plotHeatscatter(., filename = "test.pdf")