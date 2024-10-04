

# Load in libraries -------------------------------------------------------



library(here)
library(magrittr)
library(tidyverse)





# -----




# FUNCTION: Bar plot (overlay)
plotBarplotOverlay <- function(INPUT, filename){
  plotfile <- here(paste0("plots/", filename)) %T>% pdf(h=8, w=12)

  # ggplot
  rwplot <- INPUT %>%
    ggplot(data = ., aes(x = sample, y = value, fill = step)) +
    geom_bar(position = "identity", stat = "identity") +
    labs(title = "Library QC",
         subtitle = "Number of reads at each step") +
    ylab("Count") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}


# FUNCTION: Bar plot (stacked)
plotBarplotStacked <- function(INPUT, filename){
  plotfile <- here(paste0("plots/", filename)) %T>% pdf(h=8, w=12)

  # ggplot
  rwplot <- INPUT %>%
    ggplot(data = ., aes(x = sample, y = value, fill = step)) +
    geom_bar(position = "stack", stat = "identity") +
    labs(title = "Library QC",
         subtitle = "Number of reads at each step") +
    ylab("Count") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}


# FUNCTION: Bar plot (fraction)
plotBarplotFraction <- function(INPUT, filename){
  plotfile <- here(paste0("plots/", filename)) %T>% pdf(h=8, w=12)

  # ggplot
  rwplot <- INPUT %>%
    ggplot(data = ., aes(x = sample, y = value, fill = step)) +
    geom_bar(position = "fill", stat = "identity") +
    labs(title = "Library QC",
         subtitle = "Fraction of reads at each step") +
    ylab("Count") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
} 


# Read in data
readstats <- read.csv("results/post/readstats.csv", header = TRUE)


# Prepare data for plotting
INPUT <- readstats %>%
  select(-c(trim_reads)) %>% # remove, redundant with collapse_reads
  pivot_longer(cols = c(raw_reads, cutadapt_reads, collapse_reads, filter_reads, reads_mapped), names_to = "step") %>%
  mutate(step = fct_relevel(step, "raw_reads", "cutadapt_reads", "collapse_reads", "filter_reads", "reads_mapped"))
  


# Call plotting function
INPUT %>%
  plotBarplotOverlay(., filename = "barplot_readstats_overlay.pdf")


# Call plotting function
INPUT %>% mutate(step = fct_rev(step)) %>%
  plotBarplotStacked(., filename = "barplot_readstats_stacked.pdf")


# Call plotting function
INPUT %>% mutate(step = fct_rev(step)) %>%
  plotBarplotFraction(., filename = "barplot_readstats_fraction.pdf")




