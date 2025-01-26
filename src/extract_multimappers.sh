#!/bin/bash

# conda activate ribopipe

samples=("ribo_01_eIF3d_minusAux_minusHarr_4h_rep1_R1" "ribo_02_eIF3d_minusAux_minusHarr_4h_rep2_R1" "ribo_03_eIF3d_minusAux_minusHarr_4h_rep3_R1" "ribo_04_eIF3d_plusAux_minusHarr_4h_rep1_R1" "ribo_05_eIF3d_plusAux_minusHarr_4h_rep2_R1" "ribo_06_eIF3d_plusAux_minusHarr_4h_rep3_R1")

echo "Extracting multimapped reads..."

for sample in "${samples[@]}"; do
    samtools view -h "results/split_bam/transcriptome/human/${sample}.bam" | \
    awk '$1 ~ /^@/ || ($0 ~ /NH:i:[2-9]/ || $0 ~ /NH:i:[1-9][0-9]+/)' | \
    samtools view -b > "results/multimapped_bam/transcriptome/human/multimapped_${sample}.bam"
    samtools index "results/multimapped_bam/transcriptome/human/multimapped_${sample}.bam"
done

echo "Extracting unique reads..."


for sample in "${samples[@]}"; do
    samtools view -h "results/split_bam/transcriptome/human/${sample}.bam" | \
    awk '$1 ~ /^@/ || $0 ~ /NH:i:1(\s|$)/' | \
    samtools view -b > "results/unique_bam/transcriptome/human/unique_${sample}.bam"
    samtools index "results/unique_bam/transcriptome/human/unique_${sample}.bam"
done

echo "Extracting primary reads..."

for sample in "${samples[@]}"; do
    # Extract primary alignments and create the output BAM file
    samtools view -h -b -F 256 "results/split_bam/transcriptome/human/${sample}.bam" > "results/primary_bam/transcriptome/human/primary_${sample}.bam" 
    samtools index "results/primary_bam/transcriptome/human/primary_${sample}.bam"
done
