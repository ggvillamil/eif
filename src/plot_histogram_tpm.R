library(tximport)
library(ggplot2)
library(tidyverse)



gtf_gr <- readRDS("results/post/gtf_gr.rds")
txidgname <- gtf_gr %>% as.data.frame() %>% select(transcript_id, gene_name) %>% distinct()

# Pull counts (from RiboStan and Salmon)
txi_original <- readRDS("results/post/txi_conditionAuxin_human_morf_eIF3d_minusHarr_4h.rds")


tpm_df <- txi_original$abundance %>% as_tibble(rownames = "transcript_id") %>% select(-contains("ribo"))

# tpm_df <- tpm_df %>%
#   as_tibble(rownames = "transcript_id") %>%
#   select(-contains("ribo")) %>%
#   mutate(tpm_mean_minus = rowMeans(select(., contains("minusAux")))) %>%
#   mutate(tpm_mean_plus = rowMeans(select(., contains("plusAux")))) %>%
#   mutate(tpm_mean_all = rowMeans(select(., contains("total"))))



# Assuming the table has a column named "average_TPM"
plot_object <- ggplot(tpm_df, aes(x = tpm_mean_all)) +
  geom_histogram(bins = 100, fill = "steelblue", alpha = 0.7) +  # Histogram
  geom_density(color = "red", size = 1) +  # Density curve
  scale_x_log10() +  # Log scale for better visualization
  theme_minimal() +
  labs(
    title = "Distribution of Average TPM Values",
    x = "Average TPM (log10 scale)",
    y = "Density"
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black")  # Example threshold

# Save the plot
ggsave(plot_object, filename = "plots/histogram_tpm_distribution_all.pdf", width = 6, height = 4)



# Plot histogram and density curve, zooming in on TPM ≤ 10
plot_object <- ggplot(tpm_df, aes(x = tpm_mean_all)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", alpha = 0.7) +  # Finer bins for better resolution
  geom_density(color = "red", size = 1) +  # Density curve
  scale_x_continuous(trans = "pseudo_log", limits = c(0, 10)) +  # Log-like scaling but keeps zeros visible
  theme_minimal() +
  labs(
    title = "Distribution of Average TPM Values (Zoomed In: ≤10 TPM)",
    x = "Average TPM",
    y = "Density"
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black")  # Example threshold

# Save the zoomed-in plot
ggsave(plot_object, filename = "plots/histogram_tpm_distribution_zoomed_all.pdf", width = 6, height = 4)



tpm_df <- tpm_df %>%
  dplyr::rename(rna_minusAux_1 = total_01_eIF3d_minusAux_minusHarr_4h_rep1) %>%
  dplyr::rename(rna_minusAux_2 = total_02_eIF3d_minusAux_minusHarr_4h_rep2) %>%
  dplyr::rename(rna_minusAux_3 = total_03_eIF3d_minusAux_minusHarr_4h_rep3) %>%
  dplyr::rename(rna_plusAux_1 = total_04_eIF3d_plusAux_minusHarr_4h_rep1) %>%
  dplyr::rename(rna_plusAux_2 = total_05_eIF3d_plusAux_minusHarr_4h_rep2) %>%
  dplyr::rename(rna_plusAux_3 = total_06_eIF3d_plusAux_minusHarr_4h_rep3)


# Find the smallest nonzero value across all numeric columns
min_nonzero <- tpm_df %>%
  select(-transcript_id) %>%        # Exclude the character column
  as.matrix() %>%                   # Convert to a matrix for easier computation
  .[. > 0] %>%                       # Keep only nonzero values
  min(na.rm = TRUE)                  # Get the minimum nonzero value

# Replace zeros with the smallest nonzero value in numeric columns
tpm_df <- tpm_df %>%
  mutate(across(-transcript_id, ~ ifelse(. == 0, min_nonzero, .)))

# tpm_df <- tpm_df %>%
#   mutate(across(
#     -transcript_id,  # Exclude the character column
#     ~ ifelse(. == 0, mean(.[. > 0], na.rm = TRUE), .)  # Replace zero with row mean of nonzero values
#   ))



tpm_log_df <- tpm_df %>%
  mutate(rna_minusAux_1 = log(rna_minusAux_1)) %>%
  mutate(rna_minusAux_2 = log(rna_minusAux_2)) %>%
  mutate(rna_minusAux_3 = log(rna_minusAux_3)) %>%
  mutate(rna_plusAux_1 = log(rna_plusAux_1)) %>%
  mutate(rna_plusAux_2 = log(rna_plusAux_2)) %>%
  mutate(rna_plusAux_3 = log(rna_plusAux_3))

# # Replace -Inf with 0 in all numeric columns
# tpm_log_df <- tpm_log_df %>%
#   mutate(across(where(is.numeric), ~ ifelse(. == -Inf, 0, .)))

tpm_log_long <- tpm_log_df %>%
  pivot_longer(cols = -transcript_id, names_to = "library", values_to = "logTPM")



# Plot density with separate lines for each library
plot_object <- ggplot(tpm_log_long, aes(x = logTPM, color = library)) +
  geom_density(size = 1) +  # Density plot with distinct lines
  theme_minimal() +
  labs(
    title = "Density Plot of log(TPM) Values by Library",
    x = "log(TPM)",
    y = "Density"
  ) +
  theme(
    legend.title = element_blank(),  # Remove legend title
    legend.position = "right"  # Adjust legend position
  )

# Save the plot as a PDF
ggsave(plot_object, filename = "plots/density_logTPM_perLibrary_rowAverage.pdf", width = 6, height = 4)









cpm_df <- txi_original$counts



cpm_df <- cpm_df %>%
  as_tibble(rownames = "transcript_id") %>%
  select(-contains("ribo")) %>%
  dplyr::rename(rna_minusAux_1 = total_01_eIF3d_minusAux_minusHarr_4h_rep1) %>%
  dplyr::rename(rna_minusAux_2 = total_02_eIF3d_minusAux_minusHarr_4h_rep2) %>%
  dplyr::rename(rna_minusAux_3 = total_03_eIF3d_minusAux_minusHarr_4h_rep3) %>%
  dplyr::rename(rna_plusAux_1 = total_04_eIF3d_plusAux_minusHarr_4h_rep1) %>%
  dplyr::rename(rna_plusAux_2 = total_05_eIF3d_plusAux_minusHarr_4h_rep2) %>%
  dplyr::rename(rna_plusAux_3 = total_06_eIF3d_plusAux_minusHarr_4h_rep3)


total_counts <- cpm_df %>% select(-transcript_id) %>% summarise_all(.fun = sum)

cpm_df <- cpm_df %>% rowwise() %>% mutate(c_across(rna_minusAux_1:rna_plusAux_3) / total_counts, .keep = "unused") %>% ungroup()


cpm_df <- cpm_df %>%
  mutate(rna_minusAux_1 = rna_minusAux_1 * 1000000000) %>%
  mutate(rna_minusAux_2 = rna_minusAux_2 * 1000000000) %>%
  mutate(rna_minusAux_3 = rna_minusAux_3 * 1000000000) %>%
  mutate(rna_plusAux_1 = rna_plusAux_1 * 1000000000) %>%
  mutate(rna_plusAux_2 = rna_plusAux_2 * 1000000000) %>%
  mutate(rna_plusAux_3 = rna_plusAux_3 * 1000000000)


# # Find the smallest nonzero value across all numeric columns
# min_nonzero <- cpm_df %>%
#   select(-transcript_id) %>%        # Exclude the character column
#   as.matrix() %>%                   # Convert to a matrix for easier computation
#   .[. > 0] %>%                       # Keep only nonzero values
#   min(na.rm = TRUE)                  # Get the minimum nonzero value

# # Replace zeros with the smallest nonzero value in numeric columns
# cpm_df <- cpm_df %>%
#   mutate(across(-transcript_id, ~ ifelse(. == 0, min_nonzero, .)))



cpm_df <- cpm_df %>%
  mutate(across(
    -transcript_id,  # Exclude the character column
    ~ ifelse(. == 0, mean(.[. > 0], na.rm = TRUE), .)  # Replace zero with row mean of nonzero values
  ))



cpm_log_df <- cpm_df %>%
  mutate(rna_minusAux_1 = log(rna_minusAux_1)) %>%
  mutate(rna_minusAux_2 = log(rna_minusAux_2)) %>%
  mutate(rna_minusAux_3 = log(rna_minusAux_3)) %>%
  mutate(rna_plusAux_1 = log(rna_plusAux_1)) %>%
  mutate(rna_plusAux_2 = log(rna_plusAux_2)) %>%
  mutate(rna_plusAux_3 = log(rna_plusAux_3))

# # Replace -Inf with 0 in all numeric columns
# cpm_log_df <- cpm_log_df %>%
#   mutate(across(where(is.numeric), ~ ifelse(. == -Inf, 0, .)))

cpm_log_long <- cpm_log_df %>%
  pivot_longer(cols = -transcript_id, names_to = "library", values_to = "logCPM")



# Plot density with separate lines for each library
plot_object <- ggplot(cpm_log_long, aes(x = logCPM, color = library)) +
  geom_density(size = 1) +  # Density plot with distinct lines
  theme_minimal() +
  labs(
    title = "Density Plot of log(CPM) Values by Library",
    x = "log(CPM)",
    y = "Density"
  ) +
  theme(
    legend.title = element_blank(),  # Remove legend title
    legend.position = "right"  # Adjust legend position
  )

# Save the plot as a PDF
ggsave(plot_object, filename = "plots/density_logCPM_perLibrary_rowAverage.pdf", width = 6, height = 4)



