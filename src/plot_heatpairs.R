library(here)
library(magrittr)
library(tidyverse)
library(LSD) # heatscatter plot

# Get sample paths
path <- c(here("results/ribostan") %>% list.dirs(full.names = TRUE, recursive = FALSE))

# Construct sample data frame
sampledf <- tibble(path) %>%
  mutate(sample = strsplit(path, split = "/") %>% sapply(., FUN = tail, n = 1)) %>%
  separate(sample, sep = "_", into = c("assay", "subunit", "treatment", "replicate", "sample_ID", "read"), remove = FALSE) %>%
  mutate(prefix = paste0(subunit,"_",treatment))


# Functions ---------------------------------------------------------------


# FUNCTION path2sample: Extract sample name from file path
path2sample <- function(path){
  strsplit(path, split = "/") %>%
    unlist() %>%
    tail(n=1)
}


# FUNCTION read_riboem: Calculate TPMs from RiboEM (ribotrans?) output
read_riboem <- function(path){
  sample <- as.character(path) %>% strsplit(., split = "/") %>% sapply(., FUN = tail, n = 1)
  read.table(paste0(path, "/", sample, ".ribostan.human_morf_quant.tsv"), header = TRUE) %>%
    select(Name, ritpm) %>%
    rename_with(.fn = ~ paste0(sample, ".TPM"), .cols = ritpm)
}


# FUNCTION make_count_mat: Make matrix of counts (TPMs) given an input of paths
make_count_mat <- function(sampledf){
  sampledf$path %>%
  map(read_riboem) %>%
  plyr::join_all(by = "Name", type = "left")
}


# FUNCTION plot_heatpairs: Plot heatpairs given a count matrix
plot_heatpairs <- function(count_mat, main, filename, h=6, w=6){
# plot_heatpairs <- function(count_mat, h=6, w=6){
  count_mat <- count_mat %>%
    select(-Name) %>%
    mutate_all(~replace(., is.na(.), 0)) %>% # replace NAs with 0
    as.matrix %>%
    + 1 # pseudocount

  # Comment out this part to supply plot title (main) and file name (filename) manually
  {
    main <- colnames(count_mat)[1] %>%
      strsplit(split = "_") %>%
      unlist() %>%
      head(n = 3) %>% # use prefix for plot title
      # .[2] %>% # use cell-line for plot title
      paste(collapse = "_")

    filename <- paste0("heatpairs_",main,".png")
  }

  plotfile <- here(paste0("plots/",filename)) %T>%
    # pdf(h=h, w=w)
    png(h=960, w=960)

  rwplot <- heatpairs(count_mat,
    log = "xy",
    main = main)

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}


# Do the stuff ------------------------------------------------------------


# Make heatpair plots of replicates from an experiment
sampledf %>% plyr::dlply("prefix", make_count_mat) %>%
  map(plot_heatpairs)

# Compare equivalent RP and seRP experiments (compare by cell-line)
sampledf %>% filter(cellline != "helaGreen" & cellline != "helaKim") %>%
  plyr::dlply("cellline", make_count_mat) %>%
  map(plot_heatpairs, h=12, w=12)

# Compare analogous HeLa experiments
target = c("rp_hela", "rp_helaGreen", "rp_helaKim")
sampledf %>% filter(prefix %in% target) %>%
  make_count_mat %>%
  plot_heatpairs(main = "HeLa Experiments", filename = "heatpairs_helaExperiments.png", h=12, w=12)

# Compare all seRP experiments
sampledf %>% filter(assay == "serp") %>%
  make_count_mat %>%
  plot_heatpairs(main = "seRP Experiments", filename = "heatpairs_seRPExperiments.png", h=12, w=12)











##### SCRATCH




target = c("rp_hela", "rp_helaGreen", "rp_helaKim")
count_mat <- sampledf %>%
  filter(prefix %in% target) %>%
  make_count_mat

df_lists <- count_mat %>%
  tibble %>%
  select(-Name) %>%
  summarise_all(list) %>%
  pivot_longer(cols = everything(), names_to = "var", values_to = "vector")

df_lists_comb <- expand(df_lists,
  nesting(var, vector),
  nesting(var2 = var, vector2 = vector))

df_lists_comb <- df_lists_comb %>%
  filter(var != var2) %>%
  arrange(var, var2) %>%
  mutate(vars = paste0(var, ".", var2)) %>%
  select(contains("var"), everything())

c_sort_collapse <- function(...){
  c(...) %>%
    sort() %>%
    str_c(collapse = ".")
}

df_lists_comb_as <- df_lists_comb %>%
  mutate(vars = map2_chr(.x = var, 
    .y = var2,
    .f = c_sort_collapse)) %>%
  distinct(vars, .keep_all = TRUE)

pairs_cor <- df_lists_comb_as %>% 
  mutate(cor = map2(vector, vector2, cor.test, method = "pearson") %>% map_dbl("estimate"),
         vars = fct_reorder(vars, -cor))

plotfile <- here(paste0("plots/testOnly_bargraph_correlations.pdf")) %T>% pdf(h=6, w=6)

rwplot <- pairs_cor %>%
  ggplot(aes(x = vars,
    y = cor)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(title = "test",
    y = "cor",
    x = "Variable combinations") +
  theme(plot.title.position = "plot")

print(rwplot)
dev.off()
normalizePath(plotfile) %>% message


