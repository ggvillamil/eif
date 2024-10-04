
# Load libraries ----------------------------------------------------------


library(here)
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
    select(Name, NumReads)
}


# Collapse list into a data frame
readcounts_df <- bind_rows(salmon_tables, .id = "sample") %>%
  pivot_wider(names_from = sample, values_from = NumReads)


# Save read count data frame
saveRDS(readcounts_df, here("results/post/readcounts_humantotalrna.rds"))


# Build sample information table ------------------------------------------


# Specify column names
sample_cols <- c("assay", "subunit", "auxin", "replicate", "sample_id")


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
  select(-c(sample_id)) %>%
  column_to_rownames(var = "sample_name") %>%
  as.matrix()


# Save count matrix and column data tables
saveRDS(cts, here("results/post/countmatrix_humantotalrna.rds"))
saveRDS(coldata, here("results/post/deseq_coldata_totalrna.rds"))


# -----


# Optional: check if yeast transcripts have mapping reads
yeast_check <- list()
for(sample in samples){
  yeast_check[[sample]] <- read_tsv(paste0("results/salmon/data/", sample, "/quant.sf"))

  # YPL071C_mRNA is the first yeast transcript in the concatenated fasta
  yeast_check[[sample]] <- yeast_check[[sample]][which(yeast_check[[sample]]$Name == "YPL071C_mRNA"):nrow(yeast_check[[sample]]),]

  # Select non-zero read counts
  yeast_check[[sample]] <- yeast_check[[sample]][which(yeast_check[[sample]]$NumReads != 0),]
}










