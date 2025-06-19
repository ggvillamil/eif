# --- filter_transcripts.R (revised to accept args) ---

# Load libraries ----------------------------------------------------------
library(tximport)
library(tidyverse)
library(optparse)

# Argument parser ---------------------------------------------------------
option_list <- list(
  make_option("--samples", type="character", help="Space-separated list of quant.sf files", metavar="files"),
  make_option("--output", type="character", help="Output file path for transcript IDs to remove")
)

opt <- parse_args(OptionParser(option_list=option_list))

sample_files <- strsplit(opt$samples, " ")[[1]]
names(sample_files) <- basename(dirname(sample_files))

# Helper function ---------------------------------------------------------
read_quantfile <- function(filepath) {
  read_tsv(filepath) %>%
    arrange(Name) %>%
    rename_with(~ sub("ritpm", "TPM", .), .cols = contains("ritpm")) %>%
    mutate(
      TPM = replace_na(TPM, 0),
      NumReads = replace_na(NumReads, 0)
    )
}

# Import quantifications --------------------------------------------------
txi <- tximport(
  files = sample_files,
  type = "salmon",
  txOut = TRUE,
  importer = read_quantfile
)

# Load list of human transcript IDs ---------------------------------------
human_tx_bed <- read.table("resources/transcripts.human.bed", header = FALSE)
human_txid <- human_tx_bed$V1

# Filter TPM matrix -------------------------------------------------------
tpm <- txi$abundance %>% as_tibble(rownames = "transcript_id")

filtered_tpm <- tpm %>%
  filter(transcript_id %in% human_txid) %>%
  filter(rowSums(select(., -transcript_id) >= 0.01) >= 2)

remove_txid <- setdiff(human_txid, filtered_tpm$transcript_id)

# Write results -----------------------------------------------------------
write.table(
  remove_txid,
  file = opt$output,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
