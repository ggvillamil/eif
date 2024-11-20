# Load libraries ----------------------------------------------------------


library(GenomicRanges)
library(txdbmaker)
library(rtracklayer)
library(Biostrings)
library(tidyverse)


# Read in Volker's annotation ---------------------------------------------


gtf_gr <- rtracklayer::import("resources/HCT116_Txome_WT.v1.0.sort.gtf")
saveRDS(gtf_gr, "results/post/gtf_gr.rds")


# Build annotation ID table -----------------------------------------------


# Relevant columns: transcript_id, gene_id, gene_name, transcript_name, GENCODE_gene_id, GENCODE_transcript_id
anno_id_table <- mcols(gtf_gr) %>%
  as_tibble %>%
  filter(type == "transcript") %>%
  select(c(transcript_id, gene_id, gene_name, transcript_name, GENCODE_gene_id, GENCODE_transcript_id))

saveRDS(anno_id_table, "results/post/annotation_id_table.rds")
write_tsv(anno_id_table, "results/post/annotation_id_table.tsv")


# Annotate CDS start and stop ---------------------------------------------


# Make transcript database object
txdb <- makeTxDbFromGRanges(gtf_gr)
saveRDS(txdb, "results/post/gtf_txdb.rds")

# Extract exon and CDS ranges
exons_grl <- exonsBy(txdb, by = "tx", use.names = TRUE)
cds_grl <- cdsBy(txdb, by = "tx", use.names = TRUE)

# Map CDS ranges from genome space to transcript space
# Make sure the order of transcripts match between cds_grl and exons_grl
cds_trspace <- pmapToTranscripts(x = cds_grl, transcripts = exons_grl[names(cds_grl)])


# Make table of CDS start and stop positions
cds_trspace_coords <- as.data.frame(cds_trspace) %>%
  rename(transcript_id = seqnames) %>%
  mutate(transcript_id = as.character(transcript_id)) %>%
  select(c(transcript_id, start, end, width)) %>%
  rename(cds_start = start, cds_end = end, cds_width = width)

saveRDS(cds_trspace_coords, "results/post/cds_trspace_coords.rds")


# Format header string for fasta file -------------------------------------


# Relevant columns: transcript_id, transcript_name, GENCODE_transcript_id, GENCODE_transcript_name, gene_id, gene_name,  GENCODE_gene_id, GENCODE_gene_name, gene_biotype, transcript_biotype, GENCODE_gene_type, structural_category

# Extract relevant columns from annotation meta columns
anno_id_table_ext <- mcols(gtf_gr) %>%
  as_tibble %>%
  filter(type == "transcript") %>%
  select(c(transcript_id, transcript_name, GENCODE_transcript_id, GENCODE_transcript_name, gene_id, gene_name,  GENCODE_gene_id, GENCODE_gene_name, gene_biotype, transcript_biotype, GENCODE_gene_type, structural_category))

# Combine with CDS transcript space coordinates
transcript_headers_df <- left_join(anno_id_table_ext, cds_trspace_coords)

saveRDS(transcript_headers_df, "results/post/transcript_headers_df.rds")
write_tsv(transcript_headers_df, "results/post/transcript_headers_df.tsv")

# Convert to vector of strings
transcript_headers <- transcript_headers_df %>%
  unite(col = header, sep = " | ", remove = FALSE, na.rm = FALSE) %>%
  select(transcript_id, header) %>%
  deframe()


# Export transcript fasta file --------------------------------------------


# Import transcript fasta file
transcript_fasta <- rtracklayer::import("resources/HCT116_Txome_WT.v1.0.sort.fa", format = "fasta")

# Trim off CDS information from transcript fasta header names
names(transcript_fasta) <- word(names(transcript_fasta), 1)

# Replace headers with extended headers
names(transcript_fasta) <- transcript_headers[names(transcript_fasta)]

# Export as fasta file
rtracklayer::export(transcript_fasta, con = "resources/HCT116_Txome_WT.v1.0.sort.extended_header.fa", format = "fasta")


# Make CDS protein fasta file ---------------------------------------------


# Import transcript fasta file
transcript_fasta <- rtracklayer::import("resources/HCT116_Txome_WT.v1.0.sort.fa", format = "fasta")

# Trim off CDS information from transcript fasta header names
names(transcript_fasta) <- word(names(transcript_fasta), 1)

# Extract CDS sequences
cds_seqs <- Biostrings::getSeq(transcript_fasta, cds_trspace) %>%
  unlist()

# Translate to amino acid sequences
cds_aa_seqs <- Biostrings::translate(cds_seqs)

# Export as fasta file
rtracklayer::export(cds_aa_seqs, con = "resources/HCT116_Txome_WT.v1.0.sort.amino_acid.fa", format = "fasta")

# Replace headers with extended headers
names(cds_aa_seqs) <- transcript_headers[names(cds_aa_seqs)]

# Export as fasta file
rtracklayer::export(cds_aa_seqs, con = "resources/HCT116_Txome_WT.v1.0.sort.amino_acid_extended_header.fa", format = "fasta")






