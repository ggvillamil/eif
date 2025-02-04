#!/bin/bash

# Activate the required Conda environment
# Uncomment the line below if the Conda environment is not already activated
# conda activate ribopipe

# List of sample names to process
samples=(
    "total_01_eIF3d_minusAux_minusHarr_4h_rep1"
    "total_02_eIF3d_minusAux_minusHarr_4h_rep2"
    "total_03_eIF3d_minusAux_minusHarr_4h_rep3"
    "total_04_eIF3d_plusAux_minusHarr_4h_rep1"
    "total_05_eIF3d_plusAux_minusHarr_4h_rep2"
    "total_06_eIF3d_plusAux_minusHarr_4h_rep3"
)

echo "Processing reads..."

# Loop through each sample to extract multimapped reads
for sample in "${samples[@]}"; do
    # Input and output file paths
    input_bam="results/split_bam/transcriptome_rnaseq/human/${sample}.bam"
    output_bam="results/multimapped_bam/transcriptome_rnaseq/human/multimapped_${sample}.bam"

    echo "Processing ${input_file} to extract multi-mapping reads..."

    # Extract multimapped reads (NH tag > 1) and save to a new BAM file
    samtools view -h "$input_bam" | \
    awk '$1 ~ /^@/ || ($0 ~ /NH:i:[2-9]/ || $0 ~ /NH:i:[1-9][0-9]+/)' | \
    samtools view -b > "$output_bam"

    # Index the output BAM file
    samtools index "$output_bam"

    echo "Finished extracting multi-mapping reads from ${sample}. Output written to ${output_file}"
done


# Loop through each sample to extract unique reads
for sample in "${samples[@]}"; do
    # Input and output file paths
    input_bam="results/split_bam/transcriptome_rnaseq/human/${sample}.bam"
    output_bam="results/unique_bam/transcriptome_rnaseq/human/unique_${sample}.bam"

    echo "Processing ${input_file} to extract unique reads..."

    # Extract unique reads (NH tag = 1) and save to a new BAM file
    samtools view -h "$input_bam" | \
    awk '$1 ~ /^@/ || $0 ~ /NH:i:1(\s|$)/' | \
    samtools view -b > "$output_bam"

    # Index the output BAM file
    samtools index "$output_bam"

    echo "Finished extracting unique reads from ${sample}. Output written to ${output_file}"
done


# Loop through each sample to extract primary alignments
for sample in "${samples[@]}"; do
    # Input and output file paths
    input_bam="results/split_bam/transcriptome_rnaseq/human/${sample}.bam"
    output_bam="results/primary_bam/transcriptome_rnaseq/human/primary_${sample}.bam"

    echo "Processing ${input_file} to extract primary reads..."

    # Extract primary alignments (flag 256 not set) and save to a new BAM file
    samtools view -h -b -F 256 "$input_bam" > "$output_bam"

    # Index the output BAM file
    samtools index "$output_bam"

    echo "Finished extracting primary reads from ${sample}. Output written to ${output_file}"
done

echo "Processing complete."
