

# Load libraries ----------------------------------------------------------


library(here)
library(tidyverse)


# Read in RiboStan human counts -------------------------------------------


# Read sample names from sample sheet
samples <- scan(here("config/samples.csv"), skip = 1, what = "character")


# Compile read counts from RiboStan into a list with all samples
ribostan_tables <- list()
for(sample in samples){
  ribostan_tables[[sample]] <- read.table(here(paste0("results/ribostan/", sample, "/", sample, ".ribostan.human_morf_quant.tsv")), header = TRUE)

  ribostan_tables[[sample]] <- ribostan_tables[[sample]] %>%
    head(n = 43127) %>% # Remove yeast transcripts
    select(Name, NumReads) # Remove unnecessary columns
}


# Collapse list into a data frame
readcounts_df <- bind_rows(ribostan_tables, .id = "sample") %>%
  pivot_wider(names_from = sample, values_from = NumReads)


# Save read count data frame
saveRDS(readcounts_df, here("results/post/readcounts_humanribo.rds"))


# Build sample information table ------------------------------------------


# Specify column names
sample_cols <- c("assay", "subunit", "auxin", "replicate", "sample_id", "read")


# Parse column information from sample names
sample_df <- data.frame(samples) %>%
  rename(sample_name = samples) %>%
  separate(col = sample_name, into = sample_cols, sep = "_", remove = FALSE)


# Build count matrix ------------------------------------------------------


# Format read count data frame into matrix suitable for DESeq2
cts <- readcounts_df %>%
  column_to_rownames(var = "Name") %>%
  mutate_if(.predicate = is.numeric, .funs = ~ replace_na(., 0)) %>% # Replace NAs with 0
  mutate_if(.predicate = is.numeric, .funs = ~ round(., digits = 0)) %>% # Round off values to integers (RiboStan's NumReads are estimates of transcript origins)
  as.matrix()


# Format sample information table into matrix suitable for DESeq2
coldata <- sample_df %>%
  select(-c(sample_id, read)) %>%
  column_to_rownames(var = "sample_name") %>%
  as.matrix()


# Save count matrix and column data tables
saveRDS(cts, here("results/post/countmatrix_humanribo.rds"))
saveRDS(coldata, here("results/post/deseq_coldata.rds"))