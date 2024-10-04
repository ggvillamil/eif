

# Load libraries ----------------------------------------------------------


library(here)
library(magrittr)
library(DESeq2)
library(apeglm)
library(tidyverse)


# Functions ---------------------------------------------------------------


# FUNCTION: Plot MA plot
plotMAplot <- function(INPUT, filename){
  plotfile <- here(paste0("plots/", filename)) %T>% pdf(h=6, w=8)

  rwplot <- plotMA(INPUT, ylim = c(-4,4))

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}


# Select relevant samples -------------------------------------------------


# Read in count matrix, column data, and size factors
cts <- readRDS(here("results/post/countmatrix_humanribo.rds"))
coldata <- readRDS(here("results/post/deseq_coldata.rds"))
sizeFactors <- readRDS(here("results/post/sizeFactors.rds"))


# Select by subunit
subunit = "oldmonosome4G1"
i <- 25
j <- 30

my_cts <- cts[,i:j]
my_coldata <- coldata[i:j,]
my_sizeFactors <- sizeFactors[i:j]


# Build DESeq data set ----------------------------------------------------


dds <- DESeqDataSetFromMatrix(countData = my_cts,
                              colData = my_coldata,
                              design = ~ auxin)


# DESeq2 automatic size factor estimation ---------------------------------


# Run DESeq2
dds_auto <- DESeq(dds)
res_auto <- results(dds_auto)

# Bayesian shrinkage with apeglm
resLFC_auto <- lfcShrink(dds_auto, coef = "auxin_plus_vs_minus", type = "apeglm")

# Save
saveRDS(res_auto, here(paste0("results/post/deseq_res_riboOccupancy_", subunit, "_autonorm.rds")))


# Plot
res_auto %>% plotMAplot(., filename = paste0("MAplot_riboOccupancy_", subunit, "_autonorm.pdf"))
resLFC_auto %>% plotMAplot(., filename = paste0("MAplot_shrunk_riboOccupancy_", subunit, "_autonorm.pdf"))


# Size factors manually set -----------------------------------------------


# Manually set size factors
dds_manual <- dds
sizeFactors(dds_manual) <- my_sizeFactors

# Run DESeq2
dds_manual <- DESeq(dds_manual)
res_manual <- results(dds_manual)

# Bayesian shrinkage with apeglm
resLFC_manual <- lfcShrink(dds_manual, coef = "auxin_plus_vs_minus", type = "apeglm")

# Save
saveRDS(res_manual, here(paste0("results/post/deseq_res_riboOccupancy_", subunit, "_yeastnorm.rds")))

# Plot
res_manual %>% plotMAplot(., filename = paste0("MAplot_riboOccupancy_", subunit, "_yeastnorm.pdf"))
resLFC_manual %>% plotMAplot(., filename = paste0("MAplot_shrunk_riboOccupancy_", subunit, "_yeastnorm.pdf"))






