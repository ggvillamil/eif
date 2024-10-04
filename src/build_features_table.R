

# Load libraries ----------------------------------------------------------


library(here)
library(magrittr)
library(GenomicFeatures)
library(Biostrings)
library(DESeq2)
library(tidyverse)



# CCDS gtf
gtf_gr <- readRDS(here("results/post/gtf_gr.rds"))


# Calculate feature GC content and length ---------------------------------


# Make TxDb object from annotation
gtftxdb <- makeTxDbFromGRanges(gtf_gr)

# Reference genome fasta
fafile <- here("resources/GRCh38.primary_assembly.genome.fa")
fafileob <- Rsamtools::FaFile(fafile)
# Rsamtools::indexFa(fafile)

# FUNCTION: Calculate GC content and feature length
gc_content <- function(feature){
  featureseq <- extractTranscriptSeqs(fafileob, feature)

  df <- data.frame(transcript_id = names(featureseq),
    width(featureseq),
    letterFrequency(featureseq, letters = "CG", OR = "|", as.prob = TRUE))
  colnames(df)[2] = paste0("length_", substitute(feature))
  colnames(df)[3] = paste0("gc_", substitute(feature))

  df %>% as_tibble
}

# Features
cds <- cdsBy(gtftxdb, "tx", use.names = TRUE)
exons <-  exonsBy(gtftxdb, "tx", use.names = TRUE)
introns <- intronsByTranscript(gtftxdb, use.names = TRUE)
tputrs <- threeUTRsByTranscript(gtftxdb, use.names = TRUE)
fputrs <- fiveUTRsByTranscript(gtftxdb, use.names = TRUE)

startwin <- gtf_gr[gtf_gr$type =="start_codon"] %>%
  promoters(., upstream = 25, downstream = 25) %>%
  split(.$transcript_id)

stopwin <- gtf_gr[gtf_gr$type =="stop_codon"] %>%
  promoters(., upstream = 25, downstream = 25) %>%
  split(.$transcript_id)

# Apply GC content and length calculation function
gc_cds <- gc_content(cds)
gc_exons <- gc_content(exons)
gc_introns <- gc_content(introns)
gc_tputrs <- gc_content(tputrs)
gc_fputrs <- gc_content(fputrs)
gc_startwin <- gc_content(startwin)
gc_stopwin <- gc_content(stopwin)

# Fix for NaNs in intron GC content
gc_introns %<>% left_join(gc_exons) %>%
  mutate(gc_introns = if_else(is.na(gc_introns), gc_exons, gc_introns)) %>%
  dplyr::select(-c(length_exons, gc_exons))

# Join to featurestbl
featurestbl <- gc_cds %>%
  left_join(gc_exons) %>%
  left_join(gc_introns) %>%
  left_join(gc_tputrs) %>%
  left_join(gc_fputrs) %>%
  left_join(gc_startwin) %>% dplyr::select(-length_startwin) %>%
  left_join(gc_stopwin) %>% dplyr::select(-length_stopwin)

# Deal with NAs produced by joining with missing values
featurestbl %<>% replace_na(list(length_tputrs = 0, length_fputrs = 0)) %>%
  mutate(gc_tputrs = if_else(is.na(gc_tputrs), gc_exons, gc_tputrs)) %>%
  mutate(gc_fputrs = if_else(is.na(gc_fputrs), gc_exons, gc_fputrs)) %>%
  mutate(gc_startwin = if_else(is.na(gc_startwin), gc_exons, gc_startwin)) %>%
  mutate(gc_stopwin = if_else(is.na(gc_stopwin), gc_exons, gc_stopwin))








# Remove version suffix
featurestbl <- featurestbl %>%
  mutate(transcript_id_version = transcript_id) %>%
  mutate(transcript_id = sub(.$transcript_id, pattern = "\\.\\d+$", replacement = ""))


##### Motif density


densityCERT <- readRDS(here("../eif4f_archive/eif4f_2022/data/densityCERT.rds"))
densityPRTE <- readRDS(here("../eif4f_archive/eif4f_2022/data/densityPRTE.rds"))
densityTISU <- readRDS(here("../eif4f_archive/eif4f_2022/data/densityTISU.rds"))

densitySRSF1 <- readRDS(here("../eif4f_archive/eif4f_2022/data/densitySRSF1.rds"))
densityHNRPLL <- readRDS(here("../eif4f_archive/eif4f_2022/data/densityHNRPLL.rds"))
densitySAMD4A <- readRDS(here("../eif4f_archive/eif4f_2022/data/densitySAMD4A.rds"))

featurestbl <- featurestbl %>% left_join(densityCERT) %>%
  left_join(densityPRTE) %>%
  left_join(densityTISU) %>%
  left_join(densitySRSF1) %>%
  left_join(densityHNRPLL) %>%
  left_join(densitySAMD4A) %>%
  mutate(densityCERT = as.numeric(densityCERT/length_fputrs),
    densityPRTE = as.numeric(densityPRTE/length_fputrs),
    densityTISU = as.numeric(densityTISU/length_fputrs),
    densitySRSF1 = as.numeric(densitySRSF1/length_fputrs),
    densityHNRPLL = as.numeric(densityHNRPLL/length_fputrs),
    densitySAMD4A = as.numeric(densitySAMD4A/length_fputrs)) %>%
  replace_na(list(densityCERT = 0, densityPRTE = 0, densityTISU = 0, densitySRSF1 = 0, densityHNRPLL = 0, densitySAMD4A = 0))




# Save as R object
saveRDS(featurestbl, here("results/post/featurestbl.rds"))

