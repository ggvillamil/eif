library(tximport)
library(tidyverse)



gtf_gr <- readRDS("results/post/gtf_gr.rds")
txidgname <- gtf_gr %>% as.data.frame() %>% select(transcript_id, gene_name) %>% distinct()

# Pull counts (from RiboStan and Salmon)
txi_original <- readRDS("results/post/txi_conditionAuxin_human_morf_eIF3d_minusHarr_4h.rds")



tpm_df <- txi_original$abundance


te_df <- tpm_df %>%
  as_tibble(rownames = "transcript_id") %>%
  mutate(te_01_minusAux = ribo_01_eIF3d_minusAux_minusHarr_4h_rep1 / total_01_eIF3d_minusAux_minusHarr_4h_rep1) %>%
  mutate(te_02_minusAux = ribo_02_eIF3d_minusAux_minusHarr_4h_rep2 / total_02_eIF3d_minusAux_minusHarr_4h_rep2) %>%
  mutate(te_03_minusAux = ribo_03_eIF3d_minusAux_minusHarr_4h_rep3 / total_03_eIF3d_minusAux_minusHarr_4h_rep3) %>%
  mutate(te_04_plusAux = ribo_04_eIF3d_plusAux_minusHarr_4h_rep1 / total_04_eIF3d_plusAux_minusHarr_4h_rep1) %>%
  mutate(te_05_plusAux = ribo_05_eIF3d_plusAux_minusHarr_4h_rep2 / total_05_eIF3d_plusAux_minusHarr_4h_rep2) %>%
  mutate(te_06_plusAux = ribo_06_eIF3d_plusAux_minusHarr_4h_rep3 / total_06_eIF3d_plusAux_minusHarr_4h_rep3) %>%
  select(transcript_id, te_01_minusAux, te_02_minusAux, te_03_minusAux, te_04_plusAux, te_05_plusAux, te_06_plusAux)

te_log_mean_df <- te_df %>%
  mutate(te_log_mean_minusAux = log(rowMeans(select(., contains("minusAux"))))) %>%
  mutate(te_log_mean_plusAux = log(rowMeans(select(., contains("plusAux")))))






# Create scatter plot
plot_object <- ggplot(te_log_mean_df, aes(x = te_log_mean_minusAux, y = te_log_mean_plusAux)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "Scatter Plot of log(TE) Values",
    x = "log(mean(TE)) Minus Auxin",
    y = "log(mean(TE)) Plus Auxin"
  ) +
  # coord_fixed() # Ensures equal scaling for x and y axes
  xlim(-100, 100) + 
  ylim(-100, 100)

# Save the plot to a PDF file
ggsave(filename = "plots/scatter_logmeanTE_uwefilter_limits100.pdf", plot = plot_object, width = 6, height = 6)