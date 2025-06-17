
grep -vF -f results/post/table_rnaseq_filter_remove_txid_eIF3d001_4h.txt resources/HCT116_Txome_WT.v1.0.sort.gtf > resources/HCT116_Txome_WT.v1.0.sort.rnaseq_filtered_eIF3d001_4h.gtf



awk 'NR==FNR {remove[$1]; next} /^>/ {keep=!($2 in remove)} keep' results/post/table_rnaseq_filter_remove_txid_eIF3d001_4h.txt resources/HCT116_Txome_WT.v1.0.sort.fa > resources/HCT116_Txome_WT.v1.0.sort.rnaseq_filtered_eIF3d001_4h.fa

cat resources/HCT116_Txome_WT.v1.0.sort.rnaseq_filtered_eIF3d001_4h.gtf resources/Saccharomyces_cerevisiae.R64-1-1.109.gtf > resources/annotation.combined_human_yeast.rnaseq_filtered_eIF3d001_4h.gtf

cat resources/HCT116_Txome_WT.v1.0.sort.rnaseq_filtered_eIF3d001_4h.fa resources/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa > resources/transcriptome.combined_human_yeast.rnaseq_filtered_eIF3d001_4h.fa