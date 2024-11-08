

# Load libraries ----------------------------------------------------------


library(magrittr)
library(DESeq2)
library(apeglm)
library(tidyverse)


# Select samples to include in analysis -----------------------------------


my_assay <- "total" # You should replace "total" with "rnaseq" in everything
my_condition <- "conditionTime"
my_organism <- "human"
my_orf <- "morf"
my_subunit <- "eIF4G3"
my_auxin <- "plusAux"
my_harringtonine <- "minusHarr"
# my_time <- "8h"

my_suffix <- paste(my_condition, my_organism, my_orf, my_subunit, my_auxin, my_harringtonine, sep = "_")


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


sampleTable <- readRDS(paste0(paste("results/post/sampleTable", my_assay, my_suffix, sep = "_"), ".rds"))
txi <- readRDS(paste0(paste("results/post/txi", my_assay, my_suffix, sep = "_"), ".rds"))


# Build DESeq data set ----------------------------------------------------


dds <- DESeqDataSetFromTximport(txi = txi,
                              colData = sampleTable,
                              # design = ~ auxin)
                              design = ~ time)


# DESeq2 automatic size factor estimation ---------------------------------


# Run DESeq2
dds_auto <- DESeq(dds)
res_auto <- results(dds_auto)

# Save results
saveRDS(res_auto, paste0("results/post/deseq_res_diff", my_assay, "_autonorm_", my_suffix, ".rds"))

# Bayesian shrinkage with apeglm
# resLFC_auto <- lfcShrink(dds_auto, coef = "auxin_plusAux_vs_minusAux", type = "apeglm")
resLFC_auto <- lfcShrink(dds_auto, coef = "time_8h_vs_4h", type = "apeglm")

# MA plots
res_auto %>% plotMAplot(., ylim = c(-8, 8), filename = paste0("MAplot_diff", my_assay, "_autonorm_", my_suffix, ".pdf"))
resLFC_auto %>% plotMAplot(., ylim = c(-4, 4), filename = paste0("MAplot_shrunk_diff", my_assay, "_autonorm_", my_suffix, ".pdf"))


# Size factors manually set based on yeast transcripts --------------------


if(my_assay == "ribo"){
ribo_samples <- sampleTable %>% filter(assay == "ribo") %>% .$sampleName

my_sizeFactors <- readRDS("results/post/sizeFactors_yeast_ribo.rds")
names(my_sizeFactors) <- str_sub(names(my_sizeFactors), 1, -4) # You should remove the read number (i.e. _R1) from ribo-seq sample names but this will do for now
my_sizeFactors <- my_sizeFactors[ribo_samples]

# Manually set size factors
dds_manual <- dds
sizeFactors(dds_manual) <- my_sizeFactors

# Run DESeq2
dds_manual <- DESeq(dds_manual)
res_manual <- results(dds_manual)

# Save results
saveRDS(res_manual, paste0("results/post/deseq_res_diff", my_assay, "_yeastnorm_", my_suffix, ".rds"))

# Bayesian shrinkage with apeglm
# resLFC_manual <- lfcShrink(dds_manual, coef = "auxin_plusAux_vs_minusAux", type = "apeglm")
resLFC_manual <- lfcShrink(dds_manual, coef = "time_8h_vs_4h", type = "apeglm")

# MA plots
res_manual %>% plotMAplot(., ylim = c(-8, 8), filename = paste0("MAplot_diff", my_assay, "_yeastnorm_", my_suffix, ".pdf"))
resLFC_manual %>% plotMAplot(., ylim = c(-4, 4), filename = paste0("MAplot_shrunk_diff", my_assay, "_yeastnorm_", my_suffix, ".pdf"))
}


# Trying to use ribo-seq spike-ins for RNA-seq libraries ------------------


if(my_assay == "total"){
my_sizeFactors <- readRDS("results/post/sizeFactors_yeast.rds")
my_sizeFactors[82:162] <- my_sizeFactors[1:81]

total_samples <- sampleTable %>% filter(assay == "total") %>% .$sampleName
my_sizeFactors <- my_sizeFactors[total_samples]

# Manually set size factors
dds_manual <- dds
sizeFactors(dds_manual) <- my_sizeFactors

# Run DESeq2
dds_manual <- DESeq(dds_manual)
res_manual <- results(dds_manual)

# Save results
saveRDS(res_manual, paste0("results/post/deseq_res_diff", my_assay, "_yeastnorm_", my_suffix, ".rds"))

# Bayesian shrinkage with apeglm
# resLFC_manual <- lfcShrink(dds_manual, coef = "auxin_plusAux_vs_minusAux", type = "apeglm")
resLFC_manual <- lfcShrink(dds_manual, coef = "time_8h_vs_4h", type = "apeglm")

# MA plots
res_manual %>% plotMAplot(., ylim = c(-8, 8), filename = paste0("MAplot_diff", my_assay, "_yeastnorm_", my_suffix, ".pdf"))
resLFC_manual %>% plotMAplot(., ylim = c(-4, 4), filename = paste0("MAplot_shrunk_diff", my_assay, "_yeastnorm_", my_suffix, ".pdf"))
}


# Take size factors automatically estimated by DESeq2 ---------------------


if(my_assay == "total"){
  normMatrix <- assays(dds_auto)[["avgTxLength"]]
  normMatrix <- normMatrix / exp(rowMeans(log(normMatrix)))
  my_sizeFactors <- estimateSizeFactorsForMatrix(counts(dds_auto) / normMatrix)
  saveRDS(my_sizeFactors, paste0("results/post/sizeFactors_human_total_", my_suffix, ".rds"))
}


