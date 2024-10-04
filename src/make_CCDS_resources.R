

# Load R libraries --------------------------------------------------------


library(here)
library(magrittr) # pipe operator %>%
library(GenomicAlignments)
library(GenomicFeatures)
library(tidyverse) # library for data science
library(seqinr) # allows writing fasta files


# Subset annotation to CCDS entries ---------------------------------------


# Import in human annotation as gtf file
gtf <- here("resources/gencode.v37.primary_assembly.annotation.gtf")
gtf_gr <- rtracklayer::import(con = gtf, format = "gtf")


# Select CCDS entries only
gtf_gr <- subset(gtf_gr, tag == "CCDS")

# Export CCDS entries as gtf file
rtracklayer::export(gtf_gr, here("resources/gencode.v37.primary_assembly.annotation.CCDS.gtf"))

# Save
saveRDS(gtf_gr, here("results/post/gtf_gr.rds"))


# Extract transcript IDs --------------------------------------------------


# Extract transcript IDs and save as txt file
ccds_transcripts <- gtf_gr$transcript_id %>% unique()
write(ccds_transcripts, file = here("results/post/ccds_transcripts.txt"))


# Make new fasta file from selected entries only --------------------------


# Load in previously saved transcript IDs
ccds_transcripts <- scan(here("results/post/ccds_transcripts.txt"), what = "character")


# Read in fasta file of all transcripts
fafile <- here("resources/gencode.v37.pc_transcripts.fa")
fafileob <- read.fasta(file = fafile, seqtype = "DNA", as.string = TRUE, set.attributes = FALSE, forceDNAtolowe = FALSE)


# Parse fasta sequence labels to remove unneeded information
names(fafileob) <- names(fafileob) %>% str_extract("[^|]+")

# Only use transcripts that are common to CCDS entries and transcript fasta file
common <- intersect(names(fafileob), ccds_transcripts)
ccds_fafileob <- fafileob[common]


# Write new fasta file
write.fasta(sequences = ccds_fafileob, names = names(ccds_fafileob), file.out = here("resources/gencode.v37.pc_transcripts.CCDS.fa"), as.string = TRUE)