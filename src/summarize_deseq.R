# Load libraries ----------------------------------------------------------


library(DESeq2)
library(tidyverse)


# Annotation transcript ID, name, and gene information --------------------


anno_id_table <- readRDS("results/post/annotation_id_table.rds")


# summarize DESeq2 results ------------------------------------------------


subunits <- c("eIF3d", "eIF4E", "eIF4G1", "eIF4G2", "eIF4G3")
times <- c("4h", "8h")
comparison <- ("diffribo") # Choose one: "deltaTE", "diffribo", or "difftotal"
normalization <- ("yeastnorm") # Choose one: "yeastnorm" for deltaTE and diffribo; "autonorm" for difftotal


# Nested for loop
for(subunit in subunits){
  for(time in times){
    res_tbl <- readRDS(paste0("results/post/deseq_res_", comparison, "_", normalization, "_conditionAuxin_human_morf_", subunit, "_minusHarr_", time, ".rds")) %>%
      as_tibble(rownames = "transcript_id")


    res_tbl <- res_tbl %>%
      left_join(anno_id_table) %>%
      select(c(transcript_id, transcript_name, GENCODE_transcript_id, gene_id, gene_name, GENCODE_gene_id, baseMean, log2FoldChange, padj))


    write_tsv(res_tbl, paste0("results/post/deseq_table_", comparison, "_", normalization, "_conditionAuxin_human_morf_", subunit, "_minusHarr_", time, ".tsv"))


    list_up <- res_tbl %>%
      filter(log2FoldChange > 0.5 & padj < 0.05) %>%
      select(c(transcript_id, transcript_name, GENCODE_transcript_id, gene_id, gene_name, GENCODE_gene_id))

    list_down <- res_tbl %>%
      filter(log2FoldChange < -0.5 & padj < 0.05) %>%
      select(c(transcript_id, transcript_name, GENCODE_transcript_id, gene_id, gene_name, GENCODE_gene_id))


    write_tsv(list_up, paste0("results/post/list_", comparison, "_upregulated_", normalization, "_conditionAuxin_human_morf_", subunit, "_minusHarr_", time, ".tsv"))
    write_tsv(list_down, paste0("results/post/list_", comparison, "_downregulated_", normalization, "_conditionAuxin_human_morf_", subunit, "_minusHarr_", time, ".tsv"))
  }
}

