

# Load libraries ----------------------------------------------------------


library(magrittr)
library(tximport)
library(tidyverse)



# Select samples to include in analysis -----------------------------------


my_assay <- "ribo" # You should replace "total" with "rnaseq" in everything
my_condition <- "conditionTime"
my_organism <- "human"
my_orf <- "morf"
my_subunit <- "eIF3d"
my_auxin <- "plusAux"
my_harringtonine <- "minusHarr"
# my_time <- "4h"

my_suffix <- paste(my_assay, my_condition, my_organism, my_orf, my_subunit, my_auxin, my_harringtonine, sep = "_")


# Functions ---------------------------------------------------------------


# FUNCTION: Read quantification output .tsv file and rename column when necessary
read_quantfile <- function(filepath){
  read_tsv(filepath) %>%
    filter(Name %in% my_transcripts) %>% # Extract values for transcripts in your transcript list
    arrange(Name) %>%
    rename_at(vars(contains("ritpm")), list(~ sub("ritpm", "TPM", .))) %>% # tximport() looks for "TPM" column name
    mutate(TPM = replace_na(TPM, 0)) %>% # Replace NAs in TPM with 0
    mutate(NumReads = replace_na(NumReads, 0)) # Replace NAs in NumReads with 0
}


# Select transcripts ------------------------------------------------------


ribostan_transcripts <- read_tsv("results/ribostan/ribo_01_eIF3d_minusAux_minusHarr_4h_rep1_R1/ribo_01_eIF3d_minusAux_minusHarr_4h_rep1_R1.ribostan.human_morf_quant.tsv") %>% .$Name
salmon_transcripts <- read_tsv("results/salmon/data/total_01_eIF3d_minusAux_minusHarr_4h_rep1/quant.sf") %>% .$Name
my_transcripts <- intersect(ribostan_transcripts, salmon_transcripts) # Also removes yeast transcripts from Salmon output


# Construct sample information table --------------------------------------


# Specify column names
sampleCols <- c("assay", "sample_id", "subunit", "auxin", "harringtonine", "time", "replicate")


# Retrieve all sample names from sample config files
ribo_samples <- scan("config/ribo_samples.csv", skip = 1, what = "character")
ribo_samples <- str_sub(ribo_samples, 1, -4) # You should remove the read number (i.e. _R1) from ribo-seq sample names but this will do for now
rna_samples <- scan("config/total_samples.csv", skip = 1, what = "character")
samples <- c(ribo_samples, rna_samples)

# Parse column information from sample names
sampleTable <- data.frame(samples) %>%
  rename(sampleName = samples) %>%
  separate(col = sampleName, into = sampleCols, sep = "_", remove = FALSE)


# Filter samples based on specified selection -----------------------------


sampleTable <- sampleTable %>% filter(
  assay == my_assay &
  subunit == my_subunit &
  auxin == my_auxin &
  harringtonine == my_harringtonine)
  # time == my_time)

ribo_samples <- sampleTable %>% filter(assay == "ribo") %>% .$sampleName
rna_samples <- sampleTable %>% filter(assay == "total") %>% .$sampleName
samples <- sampleTable$sampleName


# Import transcript quantification ----------------------------------------


# File paths to Ribostan results
ribo_files <- paste0("results/ribostan/", ribo_samples, "_R1/", ribo_samples, "_R1.ribostan.human_morf_quant.tsv") # Again note the use of "_R1"
names(ribo_files) <- ribo_samples

# File paths to Salmon results
rna_files <- paste0("results/salmon/data/", rna_samples, "/quant.sf")
names(rna_files) <- rna_samples

# Combine Ribostan and Salmon file paths
# quantfiles <- c(ribo_files, rna_files)
quantfiles <- ribo_files

# Import abundances with tximport()
txi <- tximport(quantfiles, type = "salmon", txOut = TRUE, importer = read_quantfile)


saveRDS(sampleTable, paste0(paste("results/post/sampleTable", my_suffix, sep = "_"), ".rds"))
saveRDS(txi, paste0(paste("results/post/txi", my_suffix, sep = "_"), ".rds"))


# PCA -----

# Gene-level abundances (requires transcript to gene table)
# ribo_txi_gene <- tximport(ribo_files, type = "salmon", tx2gene = tx2genetbl, importer = read_file)

