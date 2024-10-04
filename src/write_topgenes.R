

# Load libraries ----------------------------------------------------------


library(here)
library(magrittr)
library(DESeq2)
library(tidyverse)


# Functions ---------------------------------------------------------------




# GTF GRanges copied from 2021 run (should be the same as in 2022 run)
gtf_gr <- readRDS(here("results/post/gtf_gr.rds"))


# table with gene 2 gname transformation
txid2gname <- gtf_gr %>%
  mcols() %>%
  as.data.frame() %>%
  distinct(transcript_id, gene_name) %>%
  select(transcript_id, gene_name)


# ----------


# # RiboStan results

# # Read sample names from sample sheet
# samples <- scan(here("config/samples.csv"), skip = 1, what = "character")

# # Compile read counts from RiboStan into a list with all samples
# ribostan_tables <- list()
# for(sample in samples){
#   ribostan_tables[[sample]] <- read.table(here(paste0("results/ribostan/human/", sample, "/", sample, ".human.ribostan.tsv")), header = TRUE)

#   ribostan_tables[[sample]] <- ribostan_tables[[sample]] %>%
#     head(n = 43127) %>% # Remove yeast transcripts
#     select(Name, ritpm, NumReads) # Remove unnecessary columns
# }


# ----------


# Most differentially translated


for(subunit in c("4E", "4G1", "4G2", "4G3")){
	res_dTE <- readRDS(here(paste0("results/post/deseq_res_deltaTE_", subunit, "_yeastnorm.rds")))
	res_dTE <- as.data.frame(res_dTE) %>% mutate(transcript_id = rownames(.))


	res_dTE <- res_dTE %>%
	  left_join(txid2gname) %>%
	  relocate(c(transcript_id, gene_name))


	res_dTE_down <- res_dTE %>%
	  filter(padj < 0.05) %>%
	  filter(log2FoldChange < 0) %>%
	  arrange(log2FoldChange)
	  # head(n = 500)


	res_dTE_up <- res_dTE %>%
	  filter(padj < 0.05) %>%
	  filter(log2FoldChange > 0) %>%
	  arrange(desc(log2FoldChange))
	  # head(n = 500)


	saveRDS(res_dTE_down, here(paste0("results/post/dTE_down_", subunit, ".rds")))
	saveRDS(res_dTE_up, here(paste0("results/post/dTE_up_", subunit, ".rds")))

	write.csv(res_dTE_down, file = here(paste0("results/post/dTE_down_", subunit, ".csv")), quote = FALSE, row.names = FALSE)
	write.csv(res_dTE_up, file = here(paste0("results/post/dTE_up_", subunit, ".csv")), quote = FALSE, row.names = FALSE)

	write.table(res_dTE_down$gene_name, file = here(paste0("results/post/dTE_down_", subunit, ".txt")), quote = FALSE, row.names = FALSE, col.names = FALSE)
	write.table(res_dTE_up$gene_name, file = here(paste0("results/post/dTE_up_", subunit, ".txt")), quote = FALSE, row.names = FALSE, col.names = FALSE)
}

