

# Load in libraries -------------------------------------------------------



library(here)
library(Rsamtools)
library(GenomicAlignments)
library(magrittr)
library(tidyverse)



# Functions for processing data -------------------------------------------



# FUNCTION makeGAlignments:
makeGAlignments <- function(bamFilename, organism, reference){
  bamPath <- file.path("results", "split_bam", reference, organism, bamFilename)
  bamFile <- BamFile(bamPath)

  param <- ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE), tag = "NH")

  my_GAl <- readGAlignments(bamFile, use.names = TRUE, param = param)
  return(my_GAl)
}



# FUNCTION getStats:
getStats <- function(bamFilename, reference){
  GAl_human <- makeGAlignments(bamFilename = bamFilename, organism = "human", reference = reference)
  GAl_yeast <- makeGAlignments(bamFilename = bamFilename, organism = "yeast", reference = reference)

  # Reads that align only once (Nmax == 1)
  Nmax1_human <- GAl_human[mcols(GAl_human)$NH == 1] %>% length()
  Nmax1_yeast <- GAl_yeast[mcols(GAl_yeast)$NH == 1]  %>% length()

  # All reads (including multimappers)
  reads_human <- names(GAl_human) %>% unique()
  reads_yeast <- names(GAl_yeast) %>% unique()

  # Reads (including multimappers) that map to only human or only yeast or both
  only_human <- setdiff(reads_human, reads_yeast) %>% length()
  only_yeast <- setdiff(reads_yeast, reads_human) %>% length()
  both <- intersect(reads_human, reads_yeast) %>% length()

  my_df <- data.frame("sample_name" = bamFilename,
    "Nmax1_human" = Nmax1_human,
    "Nmax1_yeast" = Nmax1_yeast,
    "only_human" = only_human,
    "only_yeast" = only_yeast,
    "both" = both)
  return(my_df)
}



# Calculate mapping statistics and build data frame -----------------------



# Pull bam file names from sample sheet
bamFilenames <- scan("config/samples.csv", what = character(), skip = 1) %>%
  paste0(., ".bam")


stats_genome <- lapply(bamFilenames, FUN = getStats, reference = "genome") %>%
  bind_rows() %>%
  rename_with(~ paste0(.x, "_genome"), !sample_name)
saveRDS(stats_genome, "results/post/stats_genome.rds")


stats_contam <- lapply(bamFilenames, FUN = getStats, reference = "contaminants") %>%
  bind_rows() %>%
  rename_with(~ paste0(.x, "_contam"), !sample_name)
saveRDS(stats_contam, "results/post/stats_contam.rds")


samplecols <- c("assay", "subunit", "depletion", "replicate", "sample_id", "mate")


statsdf <- left_join(stats_genome, stats_contam) %>%
  separate(sample_name, remove = FALSE, sep = "_", into = samplecols) %>%
  select(-c("mate"))
saveRDS(statsdf, "results/post/statsdf.rds")



# Functions for making plots ----------------------------------------------



# FUNCTION: Bar plot
plotBarplot <- function(INPUT, subtitle, filename){
  plotfile <- here(paste0("plots/", filename)) %T>% pdf(h=6, w=8)

  # ggplot
  rwplot <- INPUT %>%
    ggplot(data = ., aes(x = mapped_to, y = value)) +
    geom_bar(position = "dodge", stat = "summary", fun = "mean", aes(fill = mapped_to)) +
    geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2, position = position_dodge(0.9)) +
    facet_grid(~ depletion, switch = "x") +
    labs(title = "Count of mapping reads",
         subtitle = subtitle) +
    ylab("Count") +
    theme_bw()

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}


# FUNCTION: Bar plot
plotBarplotStacked <- function(INPUT, subtitle, filename){
  plotfile <- here(paste0("plots/", filename)) %T>% pdf(h=6, w=4)

  # ggplot
  rwplot <- INPUT %>%
    ggplot(data = ., aes(x = depletion, y = value)) +
    geom_bar(position = "fill", stat = "identity", aes(fill = mapped_to)) +
    labs(title = "Fraction of mapping reads",
         subtitle = subtitle) +
    ylab("Count") +
    theme_bw()

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}



# Let's plot! -------------------------------------------------------------



# ----- multimappers, genome

stats_genome <- readRDS("results/post/stats_genome.rds")

INPUT <- stats_genome %>%
  rename_with(~ sub("_genome", "", .x), !sample_name) %>%
  separate(sample_name, remove = FALSE, sep = "_", into = samplecols) %>%
  select(-c(mate)) %>%
  select(-c(Nmax1_human, Nmax1_yeast)) %>%
  pivot_longer(cols = c(only_human, only_yeast, both), names_to = "mapped_to")


# Call plotting function
INPUT %>% filter(subunit == "4G3") %>% select(-c(sample_name, assay, sample_id)) %>%
  plotBarplot(., subtitle = "4G3 Genome", filename = "barplot_count_multimappers_genome_4G3.pdf")



# ----- multimappers, contam

stats_contam <- readRDS("results/post/stats_contam.rds")

INPUT <- stats_contam %>%
  rename_with(~ sub("_contam", "", .x), !sample_name) %>%
  separate(sample_name, remove = FALSE, sep = "_", into = samplecols) %>%
  select(-c(mate)) %>%
  select(-c(Nmax1_human, Nmax1_yeast)) %>%
  pivot_longer(cols = c(only_human, only_yeast, both), names_to = "mapped_to")


# Call plotting function
INPUT %>% filter(subunit == "4E") %>% select(-c(sample_name, assay, sample_id)) %>%
  plotBarplot(., subtitle = "4E Contaminants", filename = "barplot_count_multimappers_contam_4E.pdf")



# ----- uniquely mapping, genome

INPUT <- stats_genome %>%
  rename_with(~ sub("_genome", "", .x), !sample_name) %>%
  separate(sample_name, remove = FALSE, sep = "_", into = samplecols) %>%
  select(-c(mate)) %>%
  select(-c(only_human, only_yeast, both)) %>%
  pivot_longer(cols = c(Nmax1_human, Nmax1_yeast), names_to = "mapped_to")


# Call plotting function
INPUT %>% filter(subunit == "4E") %>% select(-c(sample_name, assay, sample_id)) %>%
  plotBarplot(., subtitle = "4E Genome", filename = "barplot_count_uniquemappers_genome_4E.pdf")



# ----- uniquely mapping, contam

INPUT <- stats_contam %>%
  rename_with(~ sub("_contam", "", .x), !sample_name) %>%
  separate(sample_name, remove = FALSE, sep = "_", into = samplecols) %>%
  select(-c(mate)) %>%
  select(-c(only_human, only_yeast, both)) %>%
  pivot_longer(cols = c(Nmax1_human, Nmax1_yeast), names_to = "mapped_to")


# Call plotting function
INPUT %>% filter(subunit == "4G3") %>% select(-c(sample_name, assay, sample_id)) %>%
  plotBarplot(., subtitle = "4G3 Contaminants", filename = "barplot_count_uniquemappers_contam_4G3.pdf")



# ----- uniquely mapping, genome, stacked

INPUT <- stats_genome %>%
  rename_with(~ sub("_genome", "", .x), !sample_name) %>%
  separate(sample_name, remove = FALSE, sep = "_", into = samplecols) %>%
  select(-c(mate)) %>%
  select(-c(only_human, only_yeast, both)) %>%
  select(-c(sample_name, assay, sample_id)) %>%
  group_by(subunit, depletion) %>%
  summarise(Nmax1_human = mean(Nmax1_human), Nmax1_yeast = mean(Nmax1_yeast)) %>%
  ungroup() %>%
  mutate(total = Nmax1_human + Nmax1_yeast) %>%
  mutate(Nmax1_human = Nmax1_human/total, Nmax1_yeast = Nmax1_yeast/total) %>%
  select(-c(total)) %>%
  pivot_longer(cols = c(Nmax1_human, Nmax1_yeast), names_to = "mapped_to")


# Call plotting function
INPUT %>% filter(subunit == "4E") %>%
  plotBarplotStacked(., subtitle = "4E Genome", filename = "barplot_fraction_uniquemappers_genome_4E.pdf")

