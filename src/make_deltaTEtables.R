library(here)
library(DESeq2)
library(tidyverse)


# GTF GRanges copied from 2021 run (should be the same as in 2022 run)
gtf_gr <- readRDS(here("results/post/gtf_gr.rds"))


# table with gene 2 gname transformation
txid2gname <- gtf_gr %>%
  mcols() %>%
  as.data.frame() %>%
  distinct(transcript_id, gene_name) %>%
  select(transcript_id, gene_name)








# Results from DESeq deltaTE
res_files <- list.files("results/meta/deltaTE_results")
  

for(res_file in res_files){
  res <- readRDS(here(paste0("results/meta/deltaTE_results/", res_file))) %>%
    as_tibble(rownames = "transcript_id") %>%
    # drop_na(log2FoldChange) %>%
    # drop_na(padj) %>%
    select(c(transcript_id, log2FoldChange, padj))

  res <- res %>%
    left_join(txid2gname) %>%
    relocate(gene_name, .after = transcript_id)

  res <- res %>%
    arrange(log2FoldChange, padj)

  write.table(res, file = paste0("results/meta/deltaTE_tables/", sub('\\.rds$', '', res_file), "_keepNA.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
}









# Results from DESeq diffRPF
res_files <- list.files("results/meta/diffRPF_results")
  

for(res_file in res_files){
  res <- readRDS(here(paste0("results/meta/diffRPF_results/", res_file))) %>%
    as_tibble(rownames = "transcript_id") %>%
    # drop_na(log2FoldChange) %>%
    # drop_na(padj) %>%
    select(c(transcript_id, log2FoldChange, padj))

  res <- res %>%
    left_join(txid2gname) %>%
    relocate(gene_name, .after = transcript_id)

  res <- res %>%
    arrange(log2FoldChange, padj)

  write.table(res, file = paste0("results/meta/diffRPF_tables/", sub('\\.rds$', '', res_file), ".tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
}









# Results from DESeq diffExpr
res_files <- list.files("results/meta/diffExpr_results")
  

for(res_file in res_files){
  res <- readRDS(here(paste0("results/meta/diffExpr_results/", res_file))) %>%
    as_tibble(rownames = "transcript_id") %>%
    # drop_na(log2FoldChange) %>%
    # drop_na(padj) %>%
    select(c(transcript_id, log2FoldChange, padj))

  res <- res %>%
    left_join(txid2gname) %>%
    relocate(gene_name, .after = transcript_id)

  res <- res %>%
    arrange(log2FoldChange, padj)

  write.table(res, file = paste0("results/meta/diffExpr_tables/", sub('\\.rds$', '', res_file), ".tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
}
