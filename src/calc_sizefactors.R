

# Load libraries ----------------------------------------------------------


library(here)
library(tidyverse)



# Read in RiboStan yeast counts -------------------------------------------


# Read sample names from sample sheet
ribo_samples <- scan(here("config/ribo_samples.csv"), skip = 1, what = "character")

# Compile read counts from RiboStan into a list with all samples
ribostan_tables <- list()
for(sample in ribo_samples){
  ribostan_tables[[sample]] <- read.table(here(paste0("results/ribostan/", sample, "/", sample, ".ribostan.yeast_morf_quant.tsv")), header = TRUE)
  ribostan_tables[[sample]] <- ribostan_tables[[sample]] %>%
    select(Name, NumReads) # Remove unnecessary columns
}

# Collapse list into a data frame
ribo_readcounts_df <- bind_rows(ribostan_tables, .id = "sample") %>%
  pivot_wider(names_from = sample, values_from = NumReads)

yeast_transcripts <- ribo_readcounts_df$Name


# Read in Salmon yeast counts ---------------------------------------------


# Read sample names from sample sheet
rna_samples <- scan(here("config/total_samples.csv"), skip = 1, what = "character")

# Compile read counts from RiboStan into a list with all samples
salmon_tables <- list()
for(sample in rna_samples){
  salmon_tables[[sample]] <- read.table(here(paste0("results/salmon/data/", sample, "/quant.sf")), header = TRUE)
  salmon_tables[[sample]] <- salmon_tables[[sample]] %>%
    select(Name, NumReads) %>% # Remove unnecessary columns
    filter(Name %in% yeast_transcripts) # Remove human transcripts
}

# Collapse list into a data frame
rna_readcounts_df <- bind_rows(salmon_tables, .id = "sample") %>%
  pivot_wider(names_from = sample, values_from = NumReads)


# Calculate size factors --------------------------------------------------


# Join RNA-seq read counts to Ribo-seq read counts
# readcounts_df <- left_join(ribo_readcounts_df, rna_readcounts_df)
# There are zero yeast counts in the RNA-seq libraries :O

readcounts_df <- ribo_readcounts_df

# Calculate column-wise sums
sumReads <- colSums(readcounts_df[,-1], na.rm = TRUE)

# Calculate size factor relative to first library
sizeFactors <- sumReads/sumReads[1]

# Save size factors
saveRDS(sizeFactors, here("results/post/sizeFactors_yeast_ribo.rds"))