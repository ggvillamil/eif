

library(here)


# Translation initiation genes
my_genes <- c("EIF4E",
  "EIF4E2",
  "EIF4E3",
  "EIF4G1",
  "EIF4G2",
  "EIF4G3",
  "EIF2S1",
  "EIF2A",
  "EIF5",
  "EIF5B",
  "EIF3A",
  "EIF3C",
  "EIF3D",
  "EIF4A1",
  "EIF4A2",
  "EIF4EP1",
  "EIF4EP2",
  "EEF2",
  "EIF2AK3",
  "EIF2AK4")

saveRDS(my_genes, here("results/post/genelist_initiationfactors.rds"))


# NMD
my_genes <- c("NBPF8",
  "DDIT3",
  "NBPF9",
  "PPP1R15A",
  "NBPF15",
  "CEP162",
  "FAM13B",
  "NBPF12",
  "NBPF11")

saveRDS(my_genes, here("results/post/genelist_nmd.rds"))

# Membrane-specific genes
my_genes <- scan(file = here("results/post/membrane_specific_genes.txt"), what = "character")

saveRDS(my_genes, here("results/post/genelist_membranespecific.rds"))