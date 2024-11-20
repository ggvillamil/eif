library(GenomicRanges)
library(Biostrings)
library(magrittr)
# library(ggrepel) # repulsive text annotation
# library(scales) # for comma formatting numbers
# library(ggpubr) # publication-ready plots
# library(rstatix) # t_test and wilcox_test
library(tidyverse)


# Parsing functions -------------------------------------------------------


parse_ORFs_tx <- function(ORFquant_results, category){
  ORFs_tx <- ORFquant_results$ORFs_tx


  if(category != "any"){
    ORFs_tx <- ORFs_tx[ORFs_tx$ORF_category_Tx == category]
  }


  in_order <- names(ORFs_tx)
  ORFdf <- data.frame(row.names = in_order,
                      ranges(ORFs_tx)[in_order], # ORF coordinates
                      # ORFs_tx[in_order]$region, # transcript coordinates
                      gene_id = ORFs_tx[in_order]$gene_id,
                      transcript_id = ORFs_tx[in_order]$transcript_id,
                      transcript_biotype = ORFs_tx[in_order]$transcript_biotype,
                      ORF_category_Tx = ORFs_tx[in_order]$ORF_category_Tx,
                      ORF_category_Gen = ORFs_tx[in_order]$ORF_category_Gen,
                      ORFs_pM = ORFs_tx[in_order]$ORFs_pM,
                      Protein = ORFs_tx[in_order]$Protein)

  return(ORFdf)
}


make_orfdf <- function(sample, category){
  load(paste0("results/orfquant/",sample,"/", sample, "_final_ORFquant_results"))
  ORFquant_results %>% parse_ORFs_tx(category = category)
}


# Make tables from ORFquant output ----------------------------------------


ribo_samples <- scan("config/ribo_samples.csv", skip = 1, what = "character")


# Remove irrelevant samples
ribo_samples <- ribo_samples[1:60]


# Get list of parsed ORFquant outputs by using map() over all samples
orfdf_list <- ribo_samples %>%
	map(., .f = make_orfdf, category = "uORF") %>%
	setNames(ribo_samples)


# Collapse list into one data frame
my_orfdf <- orfdf_list %>%
	bind_rows(.id = "sample") %>%
	as_tibble()


# Sample information columns
my_orfdf <- my_orfdf %>%
	separate(sample, into = c("assay", "sample_id", "subunit", "auxin", "harringtonine", "time", "replicate", "read"), sep = "_", remove = FALSE)
	# mutate(group = paste(subunit, auxin, harringtonine, time, sep = "_"))


# saveRDS(my_orfdf, here("results/post/ORFquant_uorfdf.rds"))


# Plotting functions ------------------------------------------------------


# Plot number of ORFs

# FUNCTION: Bar plot
plotBarplot <- function(INPUT, subtitle, filename){
  plotfile <- paste0("plots/", filename) %T>% pdf(h=8, w=12)

  # ggplot
  rwplot <- INPUT %>%
    ggplot(data = ., aes(x = group, y = count)) +
    geom_bar(position = "dodge", stat = "summary", fun = "mean", aes(fill = auxin)) +
    geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2, position = position_dodge(0.9), aes(color = auxin)) +
    labs(title = "Number of detected uORFs",
         subtitle = subtitle) +
    ylab("Average number of detected ORFs") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}


# Number of ORFs detected per replicate
orfs_detected <- my_orfdf %>%
  mutate(group = paste(subunit, auxin, harringtonine, time, replicate, sep = "_")) %>%
  .$group %>%
  table()


orfs_detected <- orfs_detected %>%
  enframe(name = "group", value = "count") %>%
  separate(group, into = c("subunit", "auxin", "harringtonine", "time", "replicate"), sep = "_", remove = FALSE) %>%
  mutate(group = str_sub(group, 1, -6)) %>%
  mutate(group = if_else(str_detect(group, "minusAux"), paste0(str_sub(group, 1, -4), "_0h"), group)) %>%
  mutate(count = as.numeric(count))


# Call plotting function
orfs_detected %>%
  plotBarplot(., subtitle = "uORFs", filename = "barplot_ORFquant_detected_uORFs.pdf")



#####


# FUNCTION: Jitter plot
plotJitterplot <- function(INPUT, subtitle, filename){
  plotfile <- paste0("plots/", filename) %T>% pdf(h=8, w=12)

  # ggplot
  rwplot <- INPUT %>%
    ggplot(data = ., aes(x = group, y = count)) +
    geom_bar(alpha = 0.5, position = "dodge", stat = "summary", fun = "mean", aes(fill = auxin)) +
    geom_point(size = 2, position = position_jitterdodge(jitter.width = 0.05), aes(shape = auxin)) +
    # geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2, position = position_dodge(0.9), aes(color = induced)) +
    labs(title = "Number of detected uORFs",
         subtitle = subtitle) +
    ylab("Average number of detected ORFs") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}

# Call plotting function
orfs_detected %>%
  plotJitterplot(., subtitle = "uORFs", filename = "jitterplot_ORFquant_detected_uORFs.pdf")



# Plot ORFs_pM distributions as box plots ---------------------------------


# FUNCTION: Box plot
plotBoxplot <- function(INPUT, subtitle, filename){
  plotfile <- paste0("plots/", filename) %T>% pdf(h=8, w=8)

  # ggplot
  rwplot <- INPUT %>%
    ggplot(data = ., aes(x = group, y = log10(ORFs_pM))) +
    geom_boxplot(outlier.shape = 20, notch = TRUE, notchwidth = 0.8, aes(fill = auxin)) +
    # scale_y_continuous(limits = quantile(INPUT$baselineTE, c(0.1, 0.95), na.rm = TRUE)) +
    labs(title = "Change in ORF counts upon depletion",
         subtitle = subtitle) +
    ylab("log10 ORFs pM") +
    # stat_pvalue_manual(stat.test(INPUT), label = "p = {p}", xmin = "subunit", xmax = NULL) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}


# Manual statistical test table for ggplot
stat.test <- function(INPUT){
  INPUT %>%
  group_by(subunit) %>%
  t_test(baselineTE ~ induced) %>%
  adjust_pvalue() %>%
  mutate(y.position = y.position)
}


y.position = 2.25


my_orfdf <- my_orfdf %>%
  mutate(group = paste(subunit, auxin, harringtonine, time, sep = "_")) %>%
  mutate(group = if_else(str_detect(group, "minusAux"), paste0(str_sub(group, 1, -4), "_0h"), group))


# Call plotting function
my_orfdf %>%
  plotBoxplot(., subtitle = "uORFs", filename = "boxplot_ORFquant_ORFsPM_uORFs.pdf")


#####




# Plot uORF translational efficiency!


# Read in count data ------------------------------------------------------


iso_tx_countdata <- readRDS(here("../eif4f_2021/data/iso_tx_countdata.rds"))


# Remove Ribo-Seq samples from count data
rna_samples <- iso_tx_countdata$counts %>%
  colnames %>%
  str_detect("rnaseq")

iso_tx_countdata$abundance <- iso_tx_countdata$abundance[,rna_samples]
iso_tx_countdata$counts <- iso_tx_countdata$counts[,rna_samples]
iso_tx_countdata$length <- iso_tx_countdata$length[,rna_samples]




RNAtpm <- iso_tx_countdata$abundance %>%
  as_tibble(rownames = "transcript_id") %>%
  pivot_longer(!transcript_id, names_to = "sample_id", values_to = "RNAtpm")

RNAtpm <- RNAtpm %>% separate(sample_id, into = c(NA, "subunit", "induced", "replicate", NA), sep = "_", remove = FALSE)

meanRNAtpm <- RNAtpm %>%
  group_by(transcript_id, subunit, induced) %>%
  summarise(meanRNAtpm = mean(RNAtpm)) %>% ungroup()




# Remove version suffixes from transcript IDs
my_orfdf <- my_orfdf %>%
  mutate(transcript_id = sub(transcript_id, pattern = "\\.\\d+$", replacement = ""))


TEdf <- my_orfdf %>%
  select(transcript_id, subunit, induced, replicate, ORFs_pM) %>%
  left_join(meanRNAtpm) %>%
  mutate(TE = ORFs_pM / meanRNAtpm)


saveRDS(TEdf, here("data/uorf_TEdf.rds"))




# FUNCTION: Box plot
plotBoxplot2 <- function(INPUT, subtitle, filename, y.position){
  plotfile <- here(paste0("plots/", filename)) %T>% pdf(h=6, w=6)

  # ggplot
  rwplot <- INPUT %>%
    ggplot(data = ., aes(x = subunit, y = log10(TE))) +
    geom_boxplot(outlier.shape = 20, notch = TRUE, notchwidth = 0.8, aes(fill = induced)) +
    # scale_y_continuous(limits = quantile(INPUT$baselineTE, c(0.1, 0.95), na.rm = TRUE)) +
    labs(title = "Change in translational efficiency upon depletion",
         subtitle = subtitle) +
    ylab("log10 TE") +
    stat_pvalue_manual(stat.test(INPUT, y.position), label = "p = {p}", xmin = "subunit", xmax = NULL) +
    theme_bw()

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}


# Manual statistical test table for ggplot
stat.test <- function(INPUT, y.position){
  INPUT %>%
  group_by(subunit) %>%
  wilcox_test(TE ~ induced) %>%
  adjust_pvalue() %>%
  mutate(y.position = y.position)
}


# Call plotting function
TEdf %>%
  plotBoxplot2(., subtitle = "uORFs", filename = "boxplot_TE_uORFs.pdf", y.position = 3.5)




##### uORF length!



# FUNCTION: Box plot
plotBoxplot3 <- function(INPUT, subtitle, filename){
  plotfile <- here(paste0("plots/", filename)) %T>% pdf(h=6, w=6)

  # ggplot
  rwplot <- INPUT %>%
    ggplot(data = ., aes(x = subunit, y = log10(width))) +
    geom_boxplot(outlier.shape = 20, notch = TRUE, notchwidth = 0.8, aes(fill = induced)) +
    # scale_y_continuous(limits = quantile(INPUT$baselineTE, c(0.1, 0.95), na.rm = TRUE)) +
    labs(title = "Change in ORF length distribution upon depletion",
         subtitle = subtitle) +
    ylab("log10 Length") +
    # stat_pvalue_manual(stat.test(INPUT), label = "p = {p}", xmin = "subunit", xmax = NULL) +
    theme_bw()

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}

# Call plotting function
my_orfdf %>%
  plotBoxplot3(., subtitle = "uORFs", filename = "boxplot_length_uORFs.pdf")




# FUNCTION: Histograms
plotHistogram <- function(INPUT, subtitle, filename){
  plotfile <- here(paste0("plots/", filename)) %T>% pdf(h=4, w=8)

  averages <- INPUT %>% group_by(subunit, induced) %>% summarise(group_ave_width = median(width))

  # ggplot
  rwplot <- INPUT %>%
    ggplot(data = ., aes(x = log10(width), color = induced, fill = induced)) +
    geom_histogram(position = "identity", alpha = 0.5, aes(y = ..density..)) +
    facet_grid(~subunit) +
    geom_vline(data = averages, aes(xintercept = log10(group_ave_width), color = induced), linetype = "dotted") +
    geom_text_repel(data = averages, aes(x = log10(group_ave_width), y = 1.1, label = comma(round(group_ave_width)), color = induced)) +
    labs(title = "Change in ORF length distribution upon depletion",
         subtitle = subtitle) +
    xlab("log10 Length") +
    ylab("Density")
    theme_bw()

  print(rwplot)
  dev.off()
  normalizePath(plotfile) %>% message
}

# Call plotting function
my_orfdf %>%
  plotHistogram(., subtitle = "uORFs", filename = "histogram_length_uORFs.pdf")








