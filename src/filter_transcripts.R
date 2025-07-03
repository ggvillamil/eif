# Load libraries ----------------------------------------------------------

library(tximport)
library(tidyverse)

# Parse command-line arguments manually ----------------------------------
args <- commandArgs(trailingOnly = TRUE)

# Check if we have exactly 2 arguments
if (length(args) != 2) {
  stop("Usage: Rscript filter_transcripts.R <samples> <output>")
}

sample_arg <- args[1]
output_path <- args[2]

sample_files <- strsplit(sample_arg, " ")[[1]]
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
  file = output_path,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)