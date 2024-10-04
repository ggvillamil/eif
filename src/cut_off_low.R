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
    select(Name, TPM)
}


# Collapse list into a data frame
TPM_df <- bind_rows(salmon_tables, .id = "sample") %>%
  pivot_wider(names_from = sample, values_from = TPM)


# Save read count data frame
saveRDS(TPM_df, here("results/post/TPM_humantotalrna.rds"))







# -----




# # This seems way too stringent:
# greaterthan_one <- function(x){
#   x > 1
# }

# TPM_df <- TPM_df %>% filter(if_all(starts_with("total"), greaterthan_one))





TPM_df <- TPM_df %>% rowwise() %>% mutate(rowmean = mean(c_across(starts_with("total"))))

trx_pass_cutoff <- TPM_df %>%
  filter(rowmean > 1) %>%
  .$Name

saveRDS(trx_pass_cutoff, here("results/post/trx_pass_cutoff.rds"))