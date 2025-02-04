library(tximport)
library(DESeq2)
library(tidyverse)


gtf_gr <- readRDS("results/post/gtf_gr.rds")
txidgname <- gtf_gr %>% as.data.frame() %>% select(transcript_id, gene_name) %>% distinct()

# Pull counts (from RiboStan and Salmon)
txi_original <- readRDS("results/post/txi_conditionAuxin_human_morf_eIF3d_minusHarr_4h.rds")
txi_primary <- readRDS("results/post/txi_conditionAuxin_human_morf_eIF3d_minusHarr_4h_primary.rds")
txi_unique <- readRDS("results/post/txi_conditionAuxin_human_morf_eIF3d_minusHarr_4h_unique.rds")


counts_original <- txi_original$counts
counts_primary <- txi_primary$counts
counts_unique <- txi_unique$counts

# Pull DESeq results
deseq_res <- readRDS("results/post/deseq_res_deltaTE_yeastnorm_conditionAuxin_human_morf_eIF3d_minusHarr_4h.rds")

# Arrange by most extreme log2FoldChange and extract top 15
deseq_top <- deseq_res %>% as.data.frame() %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(n=15)

# Get transcript IDs of top 15
deseq_top_tx <- rownames(deseq_top)




deseq_res_primary <- readRDS("results/post/deseq_res_deltaTE_yeastnorm_conditionAuxin_human_morf_eIF3d_minusHarr_4h_primary.rds")
deseq_res_unique <- readRDS("results/post/deseq_res_deltaTE_yeastnorm_conditionAuxin_human_morf_eIF3d_minusHarr_4h_unique.rds")

deseq_top_primary <- deseq_res_primary %>% as.data.frame() %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(n=15)

deseq_top_tx_primary <- rownames(deseq_top_primary)

deseq_top_unique <- deseq_res_unique %>% as.data.frame() %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(n=15)

deseq_top_tx_unique <- rownames(deseq_top_unique)


deseq_top_tx <- c(deseq_top_tx, deseq_top_tx_primary, deseq_top_tx_unique)

compare_original <- deseq_res[deseq_top_tx,] %>% as.data.frame() %>% select(log2FoldChange) %>% cbind(counts_original[deseq_top_tx,]) %>% rownames_to_column("transcript_id")
compare_primary <- deseq_res_primary[deseq_top_tx,] %>% as.data.frame() %>% select(log2FoldChange) %>% cbind(counts_primary[deseq_top_tx,]) %>% rownames_to_column("transcript_id")
compare_unique <- deseq_res_unique[deseq_top_tx,] %>% as.data.frame() %>% select(log2FoldChange) %>% cbind(counts_unique[deseq_top_tx,]) %>% rownames_to_column("transcript_id")


compare_original <- compare_original %>% left_join(txidgname) %>% relocate(gene_name, .after = transcript_id)
compare_primary <- compare_primary %>% left_join(txidgname) %>% relocate(gene_name, .after = transcript_id)
compare_unique <- compare_unique %>% left_join(txidgname) %>% relocate(gene_name, .after = transcript_id)



write_tsv(compare_original, "results/post/table_multimap_check_top_15_original.tsv")
write_tsv(compare_primary, "results/post/table_multimap_check_top_15_primary.tsv")
write_tsv(compare_unique, "results/post/table_multimap_check_top_15_unique.tsv")




