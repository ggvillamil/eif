

# Load libraries ----------------------------------------------------------


library(here)
library(tidyverse)


# Read in RiboStan yeast counts -------------------------------------------


# Read sample names from sample sheet
samples <- scan(here("config/samples.csv"), skip = 1, what = "character")


# Compile read counts from RiboStan into a list with all samples
ribostan_tables <- list()
for(sample in samples){
  # ribostan_tables[[sample]] <- read.table(here(paste0("../eif4f_CCDS_old/results/ribostan/yeast/", sample, "/", sample, ".yeast.ribostan.tsv")), header = TRUE)
  ribostan_tables[[sample]] <- read.table(here(paste0("results/ribostan/", sample, "/", sample, ".ribostan.yeast_morf_quant.tsv")), header = TRUE)

  ribostan_tables[[sample]] <- ribostan_tables[[sample]] %>%
    # tail(n = -43127) %>% # Remove human transcripts
    select(Name, NumReads) # Remove unnecessary columns
}


# Collapse list into a data frame
readcounts_df <- bind_rows(ribostan_tables, .id = "sample") %>%
  pivot_wider(names_from = sample, values_from = NumReads)


# Save read count data frame
saveRDS(readcounts_df, here("results/post/readcounts_yeastribo.rds"))


# Calculate size factors --------------------------------------------------


# Calculate column-wise sums
sumReads <- colSums(readcounts_df[,-1], na.rm = TRUE)

# Calculate size factor relative to first library
sizeFactors <- sumReads/sumReads[1]

# Save size factors
saveRDS(sizeFactors, here("results/post/sizeFactors.rds"))