

# Load libraries ----------------------------------------------------------


library(here)
library(magrittr)
library(DESeq2)
library(apeglm)
library(tidyverse)


# Functions ---------------------------------------------------------------


# FUNCTION: Plot MA plot
plotMAplot <- function(INPUT, ylim, filename){
  plotfile <- here(paste0("plots/", filename)) %T>% pdf(h=6, w=8)

  rwplot <- plotMA(INPUT, ylim = ylim)

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}


# Select relevant samples -------------------------------------------------


# Read in count matrix, column data, and size factors
cts_ribo <- readRDS(here("results/post/countmatrix_humanribo.rds"))
coldata_ribo <- readRDS(here("results/post/deseq_coldata.rds"))
sizeFactors_ribo <- readRDS(here("results/post/sizeFactors.rds"))
cts_totalrna <- readRDS(here("results/post/countmatrix_humantotalrna.rds"))
coldata_totalrna <- readRDS(here("results/post/deseq_coldata_totalrna.rds"))


##### Try old monosome libraries of 4G1 (03.05.2024)
cts_ribo[,7] <- cts_ribo[,25]
cts_ribo[,8] <- cts_ribo[,26]
cts_ribo[,9] <- cts_ribo[,27]
cts_ribo[,10] <- cts_ribo[,28]
cts_ribo[,11] <- cts_ribo[,29]
cts_ribo[,12] <- cts_ribo[,30]
sizeFactors_ribo[7] <- sizeFactors_ribo[25]
sizeFactors_ribo[8] <- sizeFactors_ribo[26]
sizeFactors_ribo[9] <- sizeFactors_ribo[27]
sizeFactors_ribo[10] <- sizeFactors_ribo[28]
sizeFactors_ribo[11] <- sizeFactors_ribo[29]
sizeFactors_ribo[12] <- sizeFactors_ribo[30]
##### 

# Select by subunit
subunit = "4G1"
i <- 7
j <- 12

my_cts_ribo <- cts_ribo[,i:j]
my_coldata_ribo <- coldata_ribo[i:j,]
my_sizeFactors_ribo <- sizeFactors_ribo[i:j]
my_cts_totalrna <- cts_totalrna[,i:j]
my_coldata_totalrna <- coldata_totalrna[i:j,]


##### Remove bad rna-seq replicate (03.05.2024)
badrep <- "totalrna_4G1_minus_rep3_09"
badrep_index <- which(colnames(my_cts_totalrna) %in% badrep)
my_cts_totalrna <- my_cts_totalrna[,-badrep_index]
badrep_index <- which(rownames(my_coldata_totalrna) %in% badrep)
my_coldata_totalrna <- my_coldata_totalrna[-badrep_index,]
#####


common <- intersect(rownames(my_cts_ribo), rownames(my_cts_totalrna))
my_cts <- cbind(my_cts_ribo[common,], my_cts_totalrna[common,])
my_coldata <- rbind(my_coldata_ribo, my_coldata_totalrna)


# Build DESeq data set ----------------------------------------------------


# Based on deltaTE (Chothani, 2019. Current Protocols in Molecular Biology)
dds <- DESeqDataSetFromMatrix(countData = my_cts,
                              colData = my_coldata,
                              design = ~ auxin+assay+auxin:assay)


# Change reference level to RNA-seq (totalrna)
dds$assay <- relevel(dds$assay, ref = "totalrna")


# DESeq2 automatic size factor estimation ---------------------------------


# Run DESeq2
dds_auto <- DESeq(dds)
res_auto <- results(dds_auto, name = "auxinplus.assayribo")

# Bayesian shrinkage with apeglm
resLFC_auto <- lfcShrink(dds_auto, coef = "auxinplus.assayribo", type = "apeglm")

# Save
saveRDS(res_auto, here(paste0("results/post/deseq_res_deltaTE_oldmonosome", subunit, "_autonorm.rds")))


# Plot
res_auto %>% plotMAplot(., ylim = c(-8, 8), filename = paste0("MAplot_deltaTE_oldmonosome", subunit, "_autonorm.pdf"))
resLFC_auto %>% plotMAplot(., ylim = c(-4, 4), filename = paste0("MAplot_shrunk_deltaTE_oldmonosome", subunit, "_autonorm.pdf"))


# Size factors manually set -----------------------------------------------


# Estimate RNA-seq size factors independent of ribo-seq counts
dds_rnaseq <- DESeqDataSetFromMatrix(countData = my_cts_totalrna[common,],
	                                 colData = my_coldata_totalrna,
	                                 design = ~ auxin)
res_rnaseq <- DESeq(dds_rnaseq)
my_sizeFactors_totalrna <- sizeFactors(res_rnaseq)

# Manually set size factors
dds_manual <- dds

# Estimate RNA-seq size factors with ribo-seq counts # I get fewer red points (differential transcripts) this way (which makes sense)
# dds_manual <- DESeq(dds_manual)
# my_sizeFactors_totalrna <- sizeFactors(dds_manual)[7:12]/sizeFactors(dds_manual)[7:12]

sizeFactors(dds_manual) <- c(my_sizeFactors_ribo, my_sizeFactors_totalrna)

# Run DESeq2
dds_manual <- DESeq(dds_manual)
res_manual <- results(dds_manual, name = "auxinplus.assayribo")

# Bayesian shrinkage with apeglm
resLFC_manual <- lfcShrink(dds_manual, coef = "auxinplus.assayribo", type = "apeglm")

# Save
saveRDS(res_manual, here(paste0("results/post/deseq_res_deltaTE_oldmonosome", subunit, "_yeastnorm.rds")))

# Plot
res_manual %>% plotMAplot(., ylim = c(-8, 8), filename = paste0("MAplot_deltaTE_oldmonosome", subunit, "_yeastnorm.pdf"))
resLFC_manual %>% plotMAplot(., ylim = c(-4, 4), filename = paste0("MAplot_shrunk_deltaTE_oldmonosome", subunit, "_yeastnorm.pdf"))






