library(here)
library(magrittr)
library(GenomicRanges)
library(GenomicFeatures)
library(Biostrings)
library(tidyverse)


# featurestbl <- readRDS(here("results/post/featurestbl.rds"))




# Load in annotation
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

# Alternatively:
# gc_introns$gc_introns <- if_else(is.na(gc_introns$gc_introns), gc_exons$gc_exons, gc_introns$gc_introns)

# Initiate transcript features table
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




# Extract transcript_ids and gene_ids from annotation
txid2gid <- gtf_gr %>%
  mcols %>%
  as_tibble %>%
  filter(!is.na(gene_id),!is.na(transcript_id)) %>%
  distinct(transcript_id,gene_id)


featurestbl <- featurestbl %>%
  left_join(txid2gid)


# Assign transcripts to 5'TOP and ER groups -------------------------------


# Gene name to gene ID
geneName2gid <- gtf_gr %>%
    mcols %>%
    as_tibble %>%
    filter(!is.na(gene_id),!is.na(gene_name)) %>%
    distinct(gene_name,gene_id)

# Load in list of 5'TOP genes
TOPgenes <- read.csv(file = here("results/post/TOP_genes_Thoreen_2020.csv"), header = FALSE) %>% .[,1] %>%
  enframe(name = NULL, value = "gene_name") %>%
  semi_join(geneName2gid, .)
saveRDS(TOPgenes,  here("results/post/TOPgenes.rds"))

# Load in list of ER genes
ERgenes <- read.csv(file = here("results/post/uniprot_human_TMproteins_curated.tab"), sep = "\t") %>% .[,"Gene.names...primary.."] %>%
  enframe(name = NULL, value = "gene_name") %>%
  semi_join(geneName2gid, .)
saveRDS(ERgenes,  here("results/post/ERgenes.rds"))

# Prepare for joining
TOPgenes %<>% mutate(top = 1) %>%
  dplyr::select(gene_id,top)

ERgenes %<>% mutate(er = 1) %>%
  dplyr::select(gene_id,er)

# Join to featurestbl
featurestbl %<>%
  left_join(TOPgenes) %>%
  left_join(ERgenes) %>%
  replace_na(list(top = 0, er = 0))

# top = 1: transcript is transcribed from a 5'TOP gene
# er = 1: transcript is translated at the ER (i.e. comes from the uniprot list of transmembrane proteins)


# Assign genes with IRES element ------------------------------------------


# Read in IRES list from IRESbase
irestbl <- read.table(here("results/post/human_IRES_info_IRESbase20210816.txt"), header = TRUE, sep = "\t", ) %>%
  as_tibble %>%
  dplyr::rename(ires_id = IRES.ID,
    location_hg38 = Location..hg38.,
    NCBIgene_id = Gene.ID,
    gene_symbol = Gene.symbol,
    gene_synonym = Gene.synonym,
    ires_length = IRES.length)

# NCBI gene ID to Ensembl gene ID
ncbi2ensembl <- read.csv(here("results/post/biomart_ncbi2ensembl.txt")) %>%
  as_tibble %>%
  dplyr::rename(NCBIgene_id = NCBI.gene..formerly.Entrezgene..ID,
    gene_id = Gene.stable.ID.version)

# Remove entries with missing IDs
ncbi2ensembl <- ncbi2ensembl[complete.cases(ncbi2ensembl),]

# Construct IRES gene list
ires <- irestbl %>% left_join(ncbi2ensembl) %>%
  .$gene_id %>%
  as.character %>%
  tibble(gene_id = ., ires = 1)

# Save
saveRDS(ires, here("results/post/ires.rds"))

# Join to featurestbl
featurestbl %<>% left_join(ires) %>%
  replace_na(list(ires=0))






# Save

saveRDS(featurestbl, here("results/post/featurestbl.rds"))



res_file <- "results/post/deseq_res_deltaTE_4E_yeastnorm.rds"

# Results from DESeq
res <- readRDS(here(res_file)) %>%
  as_tibble(rownames = "transcript_id") %>%
  drop_na(log2FoldChange) %>%
  drop_na(padj) %>%
  dplyr::select(transcript_id, log2FoldChange)

features_dTE <- res %>%
  left_join(featurestbl)



######



# FUNCTION: Box plot
plotBoxplot <- function(INPUT, subtitle, filename, ylim){
  plotfile <- here(paste0("plots/", filename)) %T>% pdf(h=6, w=6)

  # ggplot
  rwplot <- INPUT %>%
    ggplot(data = ., aes(x = top, y = log2FoldChange)) +
    geom_boxplot(outlier.shape = NA, notch = TRUE, notchwidth = 0.8) +
    # scale_y_continuous(limits = quantile(INPUT$baselineTE, c(0.1, 0.95), na.rm = TRUE)) +
    ylim(ylim) +
    labs(title = "Change in translational efficiency upon depletion",
         subtitle = subtitle) +
    ylab("log2FC(Translational Efficiency)") +
    # stat_pvalue_manual(stat.test(INPUT), label = "p = {p}", xmin = "subunit", xmax = NULL) +
    theme_bw()

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}


# Manual statistical test table for ggplot
stat.test <- function(INPUT){
  INPUT %>%
  group_by(subunit) %>%
  t_test(baselineTE ~ induced) %>%
  adjust_pvalue() %>%
  mutate(y.position = y.position)
}


# Pare down features table
dTE <- features_dTE %>%
  dplyr::select(transcript_id, log2FoldChange, top) %>%
  mutate(top = replace(top, top == 0, "Others")) %>%
  mutate(top = replace(top, top == 1, "5TOP"))


# Call plotting function (all genes)
y.position = 2.25
dTE %>%
  plotBoxplot(., subtitle = "5'TOP vs. others", filename = "boxplot_dTE_transcripts_5primeTOP_4E.pdf", ylim = c(-5,3))











library(LSD) # heatscatter plot


# FUNCTION plot_heatpairs: Plot heatpairs given a count matrix
plot_heatpairs <- function(INPUT, main, filename){
# plot_heatpairs <- function(count_mat, h=6, w=6){
  plotfile <- here(paste0("plots/",filename)) %T>% pdf(h=6, w=6)

  rwplot <- heatpairs(INPUT,
  	method = "pearson",
    main = main)

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}



features_dTE %>% dplyr::select(length_fputrs, log2FoldChange) %>% as.matrix() %>%
  plot_heatpairs(main = "5'UTR length vs. log2FC(TE)", filename = "heatpair_lengthfputrs_4E.pdf")
