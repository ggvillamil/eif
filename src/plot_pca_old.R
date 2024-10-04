library(here)
library(magrittr)
library(tidyverse)
library(stats)


# Get sample path
path <- here("results/ribostan/human") %>% list.dirs(full.names = TRUE, recursive = FALSE)

# Construct sample data frame
sampledf <- tibble(path) %>%
  mutate(sample = strsplit(path, split = "/") %>% sapply(., FUN = tail, n = 1))


# Functions ---------------------------------------------------------------


# FUNCTION path2sample: Extract sample name from file path
path2sample <- function(path){
  strsplit(path, split = "/") %>%
    unlist() %>%
    tail(n=1)
}


# FUNCTION sampe2prefix: Extract prefix from sample name
sample2prefix <- function(sample){
  strsplit(sample, split = "_") %>%
    unlist() %>%
    head(n = 3) %>%
    tail(n = 2) %>%
    paste(collapse = "_")
}


# FUNCTION read_ribostan: Read in calculated TPMs from RiboStan
read_ribostan <- function(path){
  read.table(paste0(path,"/",path2sample(path),".human.ribostan.tsv"), header = TRUE) %>%
    head(n = 73049) %>% 
    select(Name, NumReads) %>%
    rename_with(.fn = ~ paste0(path2sample(path)), .cols = NumReads)
}


# FUNCTION make_count_mat: Make matrix of counts (TPMs) given an input of paths
make_count_mat <- function(sampledf){
  sampledf$path %>%
  map(read_ribostan) %>%
  plyr::join_all(by = "Name", type = "left")
}


#####


# Make count matrix
count_mat <- sampledf %>% make_count_mat()

# Use column "Name" as rownames then drop it
rownames(count_mat) <- count_mat[,"Name"]
count_mat <- subset(count_mat, select = -Name)

# # Transpose matrix
# count_mat <- t(count_mat)

# Center and scale input (manually)
# count_mat <- apply(count_mat, 2, function(x) (x - mean(x)) / sd(x))

# Replace NAs with 0s
count_mat[is.na(count_mat)] <- 0


sizeFactors <- readRDS(here("results/post/sizeFactors.rds"))

# Select by subunit
my_count_mat <- count_mat[,19:24]
my_sizeFactors <- sizeFactors[19:24]
my_count_mat <- my_count_mat/my_sizeFactors


# Transpose matrix
my_count_mat <- t(my_count_mat)

# Principal component analysis
# pca <- prcomp(t(count_mat), center = TRUE, scale. = TRUE)
pca <- prcomp(my_count_mat)

# Percent variation explained by principal component
pca_var <- pca$sdev^2
pca_var_per <- round(pca_var/sum(pca_var)*100, 1)

##### Plot!



# Build data frame for ggplot
pca_data <- data.frame(sample = rownames(pca$x),
  x = pca$x[,1],
  y = pca$x[,2])

# Add column for prefix (to color by)
# Add column for assay
pca_data <- pca_data %>% mutate(prefix = as.character(sample) %>% sapply(., FUN = sample2prefix)) %>%
  mutate(subunit = strsplit(prefix, split = "_") %>% sapply(., FUN = head, n = 1)) %>%
  mutate(depletion = strsplit(prefix, split = "_") %>% sapply(., FUN = tail, n = 1))

# Plot!
plotfile <- here(paste0("plots/PCAplot_NumReads_4G3_normalized.pdf")) %T>%
  pdf(h=6, w=8)
  # png(h=960, w=960)

rwplot <- ggplot(data = pca_data, aes(x = x, y = y, label = sample, color = subunit, shape = depletion)) +
  geom_point() +
  # scale_color_brewer(palette = "Set1") +
  xlab(paste0("PC1 (", pca_var_per[1], "%)")) +
  ylab(paste0("PC2 (", pca_var_per[2], "%)")) +
  theme_bw() +
  ggtitle("PCA")

print(rwplot)
dev.off()
normalizePath(plotfile) %>% message


