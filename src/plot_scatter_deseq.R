library(ggplot2)
library(dplyr)
library(ggplot2)

# Load DESeq2 result tables
deseq_res_rna <- readRDS("results/post/deseq_res_difftotal_autonorm_conditionAuxin_human_morf_eIF3d_minusHarr_4h.rds") %>%
  as_tibble(rownames = "transcript_id")
deseq_res_ribo <- readRDS("results/post/deseq_res_diffribo_yeastnorm_conditionAuxin_human_morf_eIF3d_minusHarr_4h.rds") %>%
  as_tibble(rownames = "transcript_id")


# Ensure the two tables have a common column for merging (e.g., gene ID)
merged_data <- inner_join(deseq_res_rna, deseq_res_ribo, by = "transcript_id", suffix = c("_rna", "_ribo"))

# Create scatter plot
plot_object <- ggplot(merged_data, aes(x = log2FoldChange_rna, y = log2FoldChange_ribo)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "Scatter Plot of log2FoldChange Values",
    x = "log2FoldChange (RNA)",
    y = "log2FoldChange (ribo)"
  ) +
  coord_fixed()  # Ensures equal scaling for x and y axes
  # geom_smooth(method = "lm", col = "red", se = FALSE)  # Optional linear trend line

# Save the plot to a PDF file
ggsave(filename = "plots/scatter_log2FoldChange_ribo_vs_rna_markusfilter.pdf", plot = plot_object, width = 6, height = 6)