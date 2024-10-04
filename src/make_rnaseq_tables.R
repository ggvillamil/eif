
library(here)
library(magrittr)
library(DESeq2)
library(tidyverse)









# -----


# Read sample names from sample sheet
samples <- scan(here("config/samples_totalRNA.csv"), skip = 1, what = "character")



# Compile read counts from Salmon into a list with all samples
salmon_tables <- list()
for(sample in samples){
  salmon_tables[[sample]] <- read_tsv(paste0("results/salmon/data/", sample, "/quant.sf"))

  # Remove yeast transcripts (use human transcripts only)
  salmon_tables[[sample]] <- salmon_tables[[sample]][1:which(salmon_tables[[sample]]$Name == "YPL071C_mRNA")-1,]

  # Fix transcript names
  salmon_tables[[sample]]$Name <- salmon_tables[[sample]]$Name %>%
    str_extract("[^|]+")

  # Remove unnecessary columns
  salmon_tables[[sample]] <- salmon_tables[[sample]] %>%
    select(Name, TPM)
}


# Collapse list into a data frame
tpm_df <- bind_rows(salmon_tables, .id = "sample") %>%
  pivot_wider(names_from = sample, values_from = TPM)


# Save read count data frame
saveRDS(tpm_df, here("results/post/tpm_humantotalrna.rds"))


# -----




tpm_df <- readRDS(here("results/post/tpm_humantotalrna.rds"))



subunit <- "eIF3D"
time <- "8h"
harring <- "minusHarr"

res_rnaseq <- readRDS(here(paste0("results/post/deseq_res_rnaseq_", subunit, "_", time, "_autonorm.rds")))



tpm_df_subunit <- tpm_df %>% dplyr::select(c(Name, contains(subunit))) %>% dplyr::select(c(Name, contains(time))) %>% dplyr::select(c(Name, contains(harring)))

rna_abundance_df <- res_rnaseq %>% as_tibble(rownames = "Name") %>%
  dplyr::select(c(Name, baseMean)) %>%
  left_join(tpm_df_subunit)


write_tsv(rna_abundance_df, file = here(paste0("results/post/rna_abundance_", subunit, "_", time, ".tsv")))






# -----





# Try out tximport

library(DESeq2)

files <- c("results/salmon/data/totalrna_4E_minus_rep1_01/quant.sf",
	"results/salmon/data/totalrna_4E_minus_rep2_02/quant.sf",
	"results/salmon/data/totalrna_4E_minus_rep3_03/quant.sf",
	"results/salmon/data/totalrna_4E_plus_rep1_04/quant.sf",
	"results/salmon/data/totalrna_4E_plus_rep2_05/quant.sf",
	"results/salmon/data/totalrna_4E_plus_rep3_06/quant.sf")



# dir <- system.file("extdata", package = "tximportData")
# library(readr)
# tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
# head(tx2gene)
# txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)

txi.salmon <- tximport(files, type = "salmon", txOut = TRUE)



coldata_totalrna <- readRDS(here("results/post/deseq_coldata_totalrna.rds"))
i <- 1
j <- 6
my_coldata_totalrna <- coldata_totalrna[i:j,]

colnames(txi.salmon$counts) <- rownames(my_coldata_totalrna)

# Estimate RNA-seq size factors independent of ribo-seq counts
dds_rnaseq <- DESeqDataSetFromTximport(txi = txi.salmon,
	                                 colData = my_coldata_totalrna,
	                                 design = ~ auxin)
dds_rnaseq <- DESeq(dds_rnaseq)
# my_sizeFactors_totalrna <- sizeFactors(dds_rnaseq)

res_rnaseq <- results(dds_rnaseq)
res_rnaseq %>% as_tibble(rownames = "Name")
res_rnaseq %>% as_tibble(rownames = "Name") %>% filter(Name == "ENST00000000233.10")






# -----


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



subunit <- "4E"


rna_abundance_df <- read_tsv(here(paste0("results/post/rna_abundance_", subunit, ".tsv")))

rna_abundance_df <- rna_abundance_df %>%
  rename(transcript_id = Name) %>%
  left_join(txid2gname)

most_abundant_txs <- rna_abundance_df %>%
  group_by(gene_name) %>%
  slice(which.max(baseMean)) %>%
  ungroup() %>%
  select(gene_name, transcript_id, baseMean) %>%
  arrange(gene_name)


write_tsv(most_abundant_txs, file = here(paste0("results/post/most_abundant_txs_", subunit, ".tsv")))

