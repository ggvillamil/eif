# Load necessary libraries
library(stringr)  # For string manipulation
library(tibble)
library(tidyr)
library(dplyr)
library(magrittr)
library(ggplot2)


# -----


# Define the directory path
path <- "results/star/transcriptome/data/"  # Replace with the path to your directory

# Get the contents of the directory
samples <- list.files(path)


# Functions ---------------------------------------------------------------


extract_value <- function(sample, value_label){

  # Define text file name from sample name
  file <- paste0(path, sample, "/", sample, ".transcript_Log.final.out" )

  # Read the file into a vector of lines
  lines <- readLines(file)
  
  # Search for the line containing the string
  target_line <- grep(value_label, lines, value = TRUE)

  # Extract numeric value
  value <- as.numeric(sub(".* \\|\\t([0-9]+\\.?[0-9]*).*", "\\1", target_line))

  return(value)

}


plotBoxplot <- function(INPUT, title, filename){
  plotfile <- paste0("plots/", filename) %T>% pdf(h=8, w=14)

  # ggplot
  rwplot <- INPUT %>%
    ggplot(data = ., aes(x = sample, y = percent_multi)) +
    geom_point() +
    labs(title = title) +
    ylim(c(0, 1)) +
    ylab("Percent") +
    xlab("Sample") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}


# -----


value_labels <- c("Number of reads mapped to multiple loci", "% of reads mapped to multiple loci")


value_label <- value_labels[2]
percent_multi <- sapply(samples, FUN = extract_value, value_label = value_label)

percent_multi_df <- enframe(percent_multi) %>%
  separate(name, sep = "_", into = c("assay", "sample_id", "subunit", "auxin", "harringtonine", "time", "replicate", "read")) %>%
  rename(percent_multi = value) %>%
  mutate(percent_multi = percent_multi / 100) %>%
  mutate(sample = paste(subunit, auxin, harringtonine, time, sep = "_"))


percent_multi_df %>% plotBoxplot(title = "Percent multi-mapping ribo-seq reads", filename = "boxplot_percent_multimappers_ribo.pdf")
