
# Load libraries ----------------------------------------------------------


library(here)
library(magrittr)
library(DESeq2)
library(tximport)
library(tidyverse)



# Functions ---------------------------------------------------------------


# FUNCTION: Plot PCA plot
plotPCAplot <- function(INPUT, filename, ntop){
  plotfile <- here(paste0("plots/", filename)) %T>% pdf(h=6, w=8)

  rwplot <- plotPCA(INPUT, intgroup = c("subunit", "auxin"), ntop = ntop)

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}

# Alternate:
plotPCAplot <- function(INPUT, filename, ntop){
  plotfile <- here(paste0("plots/", filename)) %T>% pdf(h=6, w=8)

  pcaData <- plotPCA(INPUT, intgroup = c("subunit", "auxin"), ntop = ntop, returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  rwplot <- ggplot(pcaData, aes(PC1, PC2, color = subunit, shape = auxin)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed()

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}


# Transcript to gene table ------------------------------------------------


gtf_gr <- readRDS(here("results/post/gtf_gr.rds"))

tx2genetbl <- gtf_gr %>%
  mcols() %>%
  as.data.frame() %>%
  filter(!is.na(gene_id), !is.na(transcript_id)) %>%
  distinct(transcript_id, gene_id) %>%
  dplyr::select(transcript_id, gene_id)

trid2gid <- tx2genetbl %>% 
  {setNames(.$gene_id, .$transcript_id)}



# Set parameters ----------------------------------------------------------


subunit <- "all"
ntop <- 100000


# Import transcript quantification ----------------------------------------


samples <- list.files("results/ribostan")
samples <- samples[-1] # Remove the folder with old results
# samples <- samples[grep(x = samples, pattern = subunit)] # Comment out to use all samples


files <- paste0("results/ribostan/", samples, "/", samples, ".ribostan.human_morf_quant.tsv")
names(files) <- samples


# Remove bad 4G1 replicate, replace with good replicate
files[grep(x = samples, pattern = "_4G1")] <- files[grep(x = samples, pattern = "oldmonosome")]
files <- files[-grep(x = samples, pattern = "oldmonosome")]
samples <- samples[-grep(x = samples, pattern = "oldmonosome")]


filereadfunc <- function(file){
  read_tsv(file) %>% mutate(ritpm = replace_na(ritpm, 0)) %>% # Replace NAs in ritpm with 0
    mutate(NumReads = replace_na(NumReads, 0)) %>% # Replace NAs in NumReads with 0
    rename(TPM = ritpm) # tximport() looks for "TPM" column name
}

# txi <- tximport(files, type = "salmon", txOut = TRUE, importer = filereadfunc)
txi <- tximport(files, type = "salmon", tx2gene = tx2genetbl, importer = filereadfunc)


# Construct sample information table --------------------------------------


coldata <- readRDS("results/post/deseq_coldata.rds")
coldata <- coldata[samples,]


# Build DESeq data set ----------------------------------------------------


# DESeq data set
dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = coldata,
                                design = ~ auxin)


# Run DESeq2
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)


plotPCAplot(INPUT = vsd, filename = paste0("PCAplot_riboseq_byGene_", subunit, "_", ntop, ".pdf"), ntop = ntop)



