library(tximport)
library(DESeq2)
library(tidyverse)


# Pull counts (from RiboStan and Salmon)
txi <- readRDS("results/post/txi_conditionAuxin_human_morf_eIF3d_minusHarr_4h.rds")
counts_df <- txi$counts
head(counts)

# Pull DESeq results
deseq_res <- readRDS("results/post/deseq_res_deltaTE_yeastnorm_conditionAuxin_human_morf_eIF3d_minusHarr_4h.rds")
deseq_res %>% as.data.frame() %>% arrange(log2FoldChange) %>% head

counts["HCT116T0006976",]
counts["HCT116T00124976",]