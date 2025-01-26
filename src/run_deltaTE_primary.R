

# Load libraries ----------------------------------------------------------


library(magrittr)
library(DESeq2)
library(apeglm)
library(tidyverse)


# Select samples to include in analysis -----------------------------------


my_condition <- "conditionAuxin"
my_organism <- "human"
my_orf <- "morf"
my_subunit <- "eIF3d"
# my_auxin <- "minusAux"
my_harringtonine <- "minusHarr"
my_time <- "4h"

my_suffix <- paste(my_condition, my_organism, my_orf, my_subunit, my_harringtonine, my_time, sep = "_")


# Functions ---------------------------------------------------------------


# FUNCTION: Plot MA plot
plotMAplot <- function(INPUT, ylim, filename){
  plotfile <- paste0("plots/", filename) %T>% pdf(h=6, w=8)

  rwplot <- plotMA(INPUT, ylim = ylim, colSig = "red")

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}


# Load sample and transcript abundance tables -----------------------------


sampleTable <- readRDS(paste0(paste("results/post/sampleTable_primary", my_suffix, sep = "_"), ".rds"))
txi <- readRDS(paste0(paste("results/post/txi_primary", my_suffix, sep = "_"), ".rds"))


# Build DESeq data set ----------------------------------------------------


# Based on deltaTE (Chothani, 2019. Current Protocols in Molecular Biology)
dds <- DESeqDataSetFromTximport(txi = txi,
                              colData = sampleTable,
                              design = ~ auxin+assay+auxin:assay)
                              # design = ~ time+assay+time:assay)


# Change reference level to RNA-seq (totalrna)
dds$assay <- relevel(dds$assay, ref = "total")


# DESeq2 automatic size factor estimation ---------------------------------


# Run DESeq2
dds_auto <- DESeq(dds)
res_auto <- results(dds_auto, name = "auxinplusAux.assayribo")
# res_auto <- results(dds_auto, name = "time4h.assayribo")

# Save results
saveRDS(res_auto, paste0("results/post/deseq_res_deltaTE_autonorm_primary_", my_suffix, ".rds"))

# Bayesian shrinkage with apeglm
resLFC_auto <- lfcShrink(dds_auto, coef = "auxinplusAux.assayribo", type = "apeglm")
# resLFC_auto <- lfcShrink(dds_auto, coef = "time4h.assayribo", type = "apeglm")

# MA plots
res_auto %>% plotMAplot(., ylim = c(-8, 8), filename = paste0("MAplot_deltaTE_autonorm_primary_", my_suffix, ".pdf"))
resLFC_auto %>% plotMAplot(., ylim = c(-4, 4), filename = paste0("MAplot_shrunk_deltaTE_autonorm_primary_", my_suffix, ".pdf"))


# Size factors manually set based on yeast transcripts --------------------


# Estimate RNA-seq size factors independent of ribo-seq counts
my_sizeFactors_rna <- readRDS(paste0("results/post/sizeFactors_human_total_", my_suffix, ".rds"))

# Retrieve ribo-seq size factors calculated based on yeast transcripts
ribo_samples <- sampleTable %>% filter(assay == "ribo") %>% .$sampleName
my_sizeFactors_ribo <- readRDS("results/post/sizeFactors_yeast_ribo.rds")
names(my_sizeFactors_ribo) <- str_sub(names(my_sizeFactors_ribo), 1, -4) # You should remove the read number (i.e. _R1) from ribo-seq sample names but this will do for now
my_sizeFactors_ribo <- my_sizeFactors_ribo[ribo_samples]

# Manually set size factors
dds_manual <- dds
sizeFactors(dds_manual) <- c(my_sizeFactors_ribo, my_sizeFactors_rna)

# Run DESeq2
dds_manual <- DESeq(dds_manual)
res_manual <- results(dds_manual, name = "auxinplusAux.assayribo")
# res_manual <- results(dds_manual, name = "time4h.assayribo")

# Save results
saveRDS(res_manual, paste0("results/post/deseq_res_deltaTE_yeastnorm_primary_", my_suffix, ".rds"))

# Bayesian shrinkage with apeglm
res_manual <- results(dds_manual, name = "auxinplusAux.assayribo", type = "apeglm")
# resLFC_manual <- lfcShrink(dds_manual, coef = "time4h.assayribo", type = "apeglm")

# MA plots
res_manual %>% plotMAplot(., ylim = c(-8, 8), filename = paste0("MAplot_deltaTE_yeastnorm_primary_", my_suffix, ".pdf"))
resLFC_manual %>% plotMAplot(., ylim = c(-4, 4), filename = paste0("MAplot_shrunk_deltaTE_yeastnorm_primary_", my_suffix, ".pdf"))

