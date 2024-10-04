

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



# Select by subunit
subunit = "4G1"
i <- 7
j <- 12

my_cts_ribo <- cts_ribo[,i:j]
my_coldata_ribo <- coldata_ribo[i:j,]
my_sizeFactors_ribo <- sizeFactors_ribo[i:j]
my_cts_totalrna <- cts_totalrna[,i:j]
my_coldata_totalrna <- coldata_totalrna[i:j,]


common <- intersect(rownames(my_cts_ribo), rownames(my_cts_totalrna))
my_cts <- cbind(my_cts_ribo[common,], my_cts_totalrna[common,])
my_coldata <- rbind(my_coldata_ribo, my_coldata_totalrna)


# Exclude problematic replicate (do duplicates)
# Potentially problematic replicate: ribo_4G1_plus_rep1_10_R1 (high level of intergenic reads)
# Try removing replicate 1 from each


# rep_names <- c("ribo_4G1_minus_rep1_07_R1",
#   "ribo_4G1_minus_rep2_08_R1",
#   "ribo_4G1_minus_rep3_09_R1",
#   "ribo_4G1_plus_rep1_10_R1",
#   "ribo_4G1_plus_rep2_11_R1",
#   "ribo_4G1_plus_rep3_12_R1")


# rep_names <- c("totalrna_4G1_minus_rep1_07",
#   "totalrna_4G1_minus_rep2_08",
#   "totalrna_4G1_minus_rep3_09",
#   "totalrna_4G1_plus_rep1_10",
#   "totalrna_4G1_plus_rep2_11",
#   "totalrna_4G1_plus_rep3_12")


# for(i in 1:6){
#   remove <- rep_names[i]

  loop_cts <- my_cts[,str_detect(colnames(my_cts), "ribo_4G1_plus_rep1_10_R1", negate = TRUE)]
  loop_coldata <- my_coldata[!rownames(my_coldata) == "ribo_4G1_plus_rep1_10_R1",]
  loop_sizeFactors_ribo <- my_sizeFactors_ribo[!names(my_sizeFactors_ribo) == "ribo_4G1_plus_rep1_10_R1"]

  loop_cts <- loop_cts[,str_detect(colnames(loop_cts), "totalrna_4G1_minus_rep3_09", negate = TRUE)]
  loop_coldata <- loop_coldata[!rownames(loop_coldata) == "totalrna_4G1_minus_rep3_09",]
  loop_sizeFactors_ribo <- loop_sizeFactors_ribo[!names(loop_sizeFactors_ribo) == "totalrna_4G1_minus_rep3_09"]

  # my_cts <- my_cts[,str_detect(colnames(my_cts), "rep1", negate = TRUE)]
  # my_coldata <- my_coldata[my_coldata[,"replicate"] != "rep1",]
  # my_sizeFactors_ribo <- my_sizeFactors_ribo[str_detect(names(my_sizeFactors_ribo), "rep1", negate = TRUE)]


  loop_cts_totalrna <- my_cts_totalrna[,str_detect(colnames(my_cts_totalrna), "totalrna_4G1_minus_rep3_09", negate = TRUE)]
  loop_coldata_totalrna <- my_coldata_totalrna[!rownames(my_coldata_totalrna) == "totalrna_4G1_minus_rep3_09",]

  # my_cts_totalrna <- my_cts_totalrna[,str_detect(colnames(my_cts_totalrna), "rep1", negate = TRUE)]
  # my_coldata_totalrna <- my_coldata_totalrna[my_coldata_totalrna[,"replicate"] != "rep1",]


  # Build DESeq data set ----------------------------------------------------


  # Based on deltaTE (Chothani, 2019. Current Protocols in Molecular Biology)
  dds <- DESeqDataSetFromMatrix(countData = loop_cts,
                                colData = loop_coldata,
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
  saveRDS(res_auto, here(paste0("results/post/deseq_res_deltaTE_", subunit, "_autonorm_4G1check_test4-3.rds")))


  # Plot
  res_auto %>% plotMAplot(., ylim = c(-8, 8), filename = paste0("MAplot_deltaTE_", subunit, "_autonorm_4G1check_test4-3.pdf"))
  resLFC_auto %>% plotMAplot(., ylim = c(-4, 4), filename = paste0("MAplot_shrunk_deltaTE_", subunit, "_autonorm_4G1check_test4-3.pdf"))


  # Size factors manually set -----------------------------------------------


  # Estimate RNA-seq size factors independent of ribo-seq counts
  dds_rnaseq <- DESeqDataSetFromMatrix(countData = loop_cts_totalrna[common,],
                                     colData = loop_coldata_totalrna,
                                     design = ~ auxin)
  res_rnaseq <- DESeq(dds_rnaseq)
  my_sizeFactors_totalrna <- sizeFactors(res_rnaseq)

  # Manually set size factors
  dds_manual <- dds

  # Estimate RNA-seq size factors with ribo-seq counts # I get fewer red points (differential transcripts) this way (which makes sense)
  # dds_manual <- DESeq(dds_manual)
  # my_sizeFactors_totalrna <- sizeFactors(dds_manual)[7:12]/sizeFactors(dds_manual)[7:12]

  sizeFactors(dds_manual) <- c(loop_sizeFactors_ribo, my_sizeFactors_totalrna)

  # Run DESeq2
  dds_manual <- DESeq(dds_manual)
  res_manual <- results(dds_manual, name = "auxinplus.assayribo")

  # Bayesian shrinkage with apeglm
  resLFC_manual <- lfcShrink(dds_manual, coef = "auxinplus.assayribo", type = "apeglm")

  # Save
  saveRDS(res_manual, here(paste0("results/post/deseq_res_deltaTE_", subunit, "_yeastnorm_4G1check_test4-3.rds")))

  # Plot
  res_manual %>% plotMAplot(., ylim = c(-8, 8), filename = paste0("MAplot_deltaTE_", subunit, "_yeastnorm_4G1check_test4-3.pdf"))
  resLFC_manual %>% plotMAplot(., ylim = c(-4, 4), filename = paste0("MAplot_shrunk_deltaTE_", subunit, "_yeastnorm_4G1check_test4-3.pdf"))


# }



