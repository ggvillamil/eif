

# Load libraries ----------------------------------------------------------


library(tximport)
library(tidyverse)


# Functions ---------------------------------------------------------------


# FUNCTION: Read quantification output .tsv file and rename column when necessary
read_quantfile <- function(filepath){
  read_tsv(filepath) %>%
    # filter(Name %in% my_transcripts) %>% # Extract values for transcripts in your transcript list
    arrange(Name) %>%
    rename_at(vars(contains("ritpm")), list(~ sub("ritpm", "TPM", .))) %>% # tximport() looks for "TPM" column name
    mutate(TPM = replace_na(TPM, 0)) %>% # Replace NAs in TPM with 0
    mutate(NumReads = replace_na(NumReads, 0)) # Replace NAs in NumReads with 0
}


# Import transcript quantification ----------------------------------------


# Retrieve all sample names from sample config files
rna_samples <- scan("config/total_samples.csv", skip = 1, what = "character")

rna_samples <- rna_samples[25:30]


# File paths to Salmon results
rna_files <- paste0("results/salmon/data/", rna_samples, "/quant.sf")
names(rna_files) <- rna_samples

# Import abundances with tximport()
txi <- tximport(rna_files, type = "salmon", txOut = TRUE, importer = read_quantfile)


# Human transcript IDs -----


human_tx_bed <- read.table("resources/transcripts.human.bed", header = FALSE)
human_txid <- human_tx_bed$V1


# Filter -----


# Put TPM values in a tibble
tpm <- txi$abundance %>% as_tibble(rownames = "transcript_id")

# Uwe Filter: rows must have tpm >= 1 in at least 1 column in both conditions
# Remove yeast transcripts along the way
# filtered_tpm <- tpm %>%
#   filter(transcript_id %in% human_txid) %>%
#   filter(
#     rowSums(select(., contains("minusAux")) >= 1) >= 1 &
#     rowSums(select(., contains("plusAux")) >= 1) >= 1
#   )

# Markus Filter: rows must have tpm >= 1 in at least 2 columns in either condition
# Remove yeast transcripts along the way
filtered_tpm <- tpm %>%
  filter(transcript_id %in% human_txid) %>%
  filter(rowSums(select(., -transcript_id) >= 1) >= 2)


filtered_txid <- filtered_tpm$transcript_id

remove_txid <- setdiff(human_txid, filtered_txid)
# write.table(remove_txid, file = "results/post/table_rnaseq_filter_remove_txid_eIF3d_uwefilter.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(remove_txid, file = "results/post/table_rnaseq_filter_remove_txid_eIF4G3_4h.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


