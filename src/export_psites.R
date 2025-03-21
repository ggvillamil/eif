library(here)
library(GenomicRanges)
library(GenomicFeatures)
library(txdbmaker)
library(rtracklayer)
library(tidyverse)



# Scan directory for sample names
samples <- list.files("results/ribostan")


for(sample in samples){

  # Load in p-sites
  psites <- readRDS(paste0("results/ribostan/", sample, "/", sample, ".ribostan.human_psites.rds"))

  # Build transcript database
  gtf_gr <- readRDS(here("results/post/gtf_gr.rds"))
  gtftxdb <- makeTxDbFromGRanges(gtf_gr)

  # Convert transcript database to GenomicRanges (again)
  gtf_tx_gr <- exonsBy(gtftxdb, use.names = TRUE)

  # Map from transcripts to genome
  psites_genome <- mapFromTranscripts(x = psites, transcripts = gtf_tx_gr, ignore.strand = TRUE)

  # Convert to converage (RleList)
  psites_cov <- coverage(psites_genome)

  # Export as bedGraph
  rtracklayer::export(psites_cov, con = paste0("results/post/psitecov_", sample, ".bedgraph"), format = "bedGraph")

}


# Everything looks spliced! [FIXED by using exonsBy() instead of transcripts()]
# And way too sparse
# A lot of p-sites lost when converting from psites to psites_genome [FIXED by setting flag ignore.strand = TRUE in mapFromTranscripts()]
