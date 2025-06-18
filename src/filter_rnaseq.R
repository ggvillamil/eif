# -------------------------------------------------------------------------
# Filter transcript quantification based on expression thresholds
# -------------------------------------------------------------------------
# This script reads TPM data from Salmon quantification files, filters
# transcripts based on expression across samples, and outputs a list of
# transcripts to be removed from further analysis.
# -------------------------------------------------------------------------

# Load libraries ----------------------------------------------------------

library(tximport)
library(tidyverse)

# Define helper functions -------------------------------------------------

# Read a Salmon quantification .tsv file and standardize column names
read_quantfile <- function(filepath) {
  read_tsv(filepath) %>%
    arrange(Name) %>%
    rename_with(~ sub("ritpm", "TPM", .), .cols = contains("ritpm")) %>%
    mutate(
      TPM = replace_na(TPM, 0),
      NumReads = replace_na(NumReads, 0)
    )
}

# Load RNA-seq sample information -----------------------------------------

# Read sample IDs from CSV (skip header row)
rna_samples <- scan("config/total_samples.csv", skip = 1, what = "character")

# Find way to select per subunit and treatment

# Generate file paths to Salmon output
rna_files <- paste0("results/salmon/data/", rna_samples, "/quant.sf")
names(rna_files) <- rna_samples

# Import transcript-level abundance estimates -----------------------------

txi <- tximport(
  files = rna_files,
  type = "salmon",
  txOut = TRUE,
  importer = read_quantfile
)

# Load list of human transcript IDs ---------------------------------------

human_tx_bed <- read.table("resources/transcripts.human.bed", header = FALSE)
human_txid <- human_tx_bed$V1

# Prepare TPM matrix for filtering ----------------------------------------

tpm <- txi$abundance %>%
  as_tibble(rownames = "transcript_id")

# Filter transcripts ------------------------------------------------------

# Keep transcripts with TPM >= 0.01 in at least 2 samples (any condition)
filtered_tpm <- tpm %>%
  filter(transcript_id %in% human_txid) %>%
  filter(rowSums(select(., -transcript_id) >= 0.01) >= 2)

# Identify transcripts to remove ------------------------------------------

filtered_txid <- filtered_tpm$transcript_id
remove_txid <- setdiff(human_txid, filtered_txid)

# Write output ------------------------------------------------------------

write.table(
  remove_txid,
  file = "results/post/table_rnaseq_filter_remove_txid_eIF3d001_4h.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
