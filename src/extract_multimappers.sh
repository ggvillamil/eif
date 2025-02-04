#!/bin/bash

# Activate the required Conda environment
# Uncomment the line below if the Conda environment is not already activated
# conda activate ribopipe

# List of sample names to process
samples=(
    "ribo_07_eIF4E_minusAux_minusHarr_4h_rep1_R1"
    "ribo_08_eIF4E_minusAux_minusHarr_4h_rep2_R1"
    "ribo_09_eIF4E_minusAux_minusHarr_4h_rep3_R1"
    "ribo_10_eIF4E_plusAux_minusHarr_4h_rep1_R1"
    "ribo_11_eIF4E_plusAux_minusHarr_4h_rep2_R1"
    "ribo_12_eIF4E_plusAux_minusHarr_4h_rep3_R1"
    "ribo_13_eIF4G1_minusAux_minusHarr_4h_rep1_R1"
    "ribo_14_eIF4G1_minusAux_minusHarr_4h_rep2_R1"
    "ribo_15_eIF4G1_minusAux_minusHarr_4h_rep3_R1"
    "ribo_16_eIF4G1_plusAux_minusHarr_4h_rep1_R1"
    "ribo_17_eIF4G1_plusAux_minusHarr_4h_rep2_R1"
    "ribo_18_eIF4G1_plusAux_minusHarr_4h_rep3_R1"
    "ribo_19_eIF4G2_minusAux_minusHarr_4h_rep1_R1"
    "ribo_20_eIF4G2_minusAux_minusHarr_4h_rep2_R1"
    "ribo_21_eIF4G2_minusAux_minusHarr_4h_rep3_R1"
    "ribo_22_eIF4G2_plusAux_minusHarr_4h_rep1_R1"
    "ribo_23_eIF4G2_plusAux_minusHarr_4h_rep2_R1"
    "ribo_24_eIF4G2_plusAux_minusHarr_4h_rep3_R1"
    "ribo_25_eIF4G3_minusAux_minusHarr_4h_rep1_R1"
    "ribo_26_eIF4G3_minusAux_minusHarr_4h_rep2_R1"
    "ribo_27_eIF4G3_minusAux_minusHarr_4h_rep3_R1"
    "ribo_28_eIF4G3_plusAux_minusHarr_4h_rep1_R1"
    "ribo_29_eIF4G3_plusAux_minusHarr_4h_rep2_R1"
    "ribo_30_eIF4G3_plusAux_minusHarr_4h_rep3_R1"
    "ribo_31_eIF3d_minusAux_minusHarr_8h_rep1_R1"
    "ribo_32_eIF3d_minusAux_minusHarr_8h_rep2_R1"
    "ribo_33_eIF3d_minusAux_minusHarr_8h_rep3_R1"
    "ribo_34_eIF3d_plusAux_minusHarr_8h_rep1_R1"
    "ribo_35_eIF3d_plusAux_minusHarr_8h_rep2_R1"
    "ribo_36_eIF3d_plusAux_minusHarr_8h_rep3_R1"
    "ribo_37_eIF4E_minusAux_minusHarr_8h_rep1_R1"
    "ribo_38_eIF4E_minusAux_minusHarr_8h_rep2_R1"
    "ribo_39_eIF4E_minusAux_minusHarr_8h_rep3_R1"
    "ribo_40_eIF4E_plusAux_minusHarr_8h_rep1_R1"
    "ribo_41_eIF4E_plusAux_minusHarr_8h_rep2_R1"
    "ribo_42_eIF4E_plusAux_minusHarr_8h_rep3_R1"
    "ribo_43_eIF4G1_minusAux_minusHarr_8h_rep1_R1"
    "ribo_44_eIF4G1_minusAux_minusHarr_8h_rep2_R1"
    "ribo_45_eIF4G1_minusAux_minusHarr_8h_rep3_R1"
    "ribo_46_eIF4G1_plusAux_minusHarr_8h_rep1_R1"
    "ribo_47_eIF4G1_plusAux_minusHarr_8h_rep2_R1"
    "ribo_48_eIF4G1_plusAux_minusHarr_8h_rep3_R1"
    "ribo_49_eIF4G2_minusAux_minusHarr_8h_rep1_R1"
    "ribo_50_eIF4G2_minusAux_minusHarr_8h_rep2_R1"
    "ribo_51_eIF4G2_minusAux_minusHarr_8h_rep3_R1"
    "ribo_52_eIF4G2_plusAux_minusHarr_8h_rep1_R1"
    "ribo_53_eIF4G2_plusAux_minusHarr_8h_rep2_R1"
    "ribo_54_eIF4G2_plusAux_minusHarr_8h_rep3_R1"
    "ribo_55_eIF4G3_minusAux_minusHarr_8h_rep1_R1"
    "ribo_56_eIF4G3_minusAux_minusHarr_8h_rep2_R1"
    "ribo_57_eIF4G3_minusAux_minusHarr_8h_rep3_R1"
    "ribo_58_eIF4G3_plusAux_minusHarr_8h_rep1_R1"
    "ribo_59_eIF4G3_plusAux_minusHarr_8h_rep2_R1"
    "ribo_60_eIF4G3_plusAux_minusHarr_8h_rep3_R1"
    "ribo_61_HCT116_minusAux_minusHarr_8h_rep1_R1"
    "ribo_62_HCT116_minusAux_minusHarr_8h_rep2_R1"
    "ribo_63_HCT116_minusAux_minusHarr_8h_rep3_R1"
    "ribo_64_HCT116_plusAux_minusHarr_8h_rep1_R1"
    "ribo_65_HCT116_plusAux_minusHarr_8h_rep2_R1"
    "ribo_66_HCT116_plusAux_minusHarr_8h_rep3_R1"
    "ribo_67_eIF3d_plusAux_plusHarr_4h_rep1_R1"
    "ribo_68_eIF3d_plusAux_plusHarr_4h_rep2_R1"
    "ribo_69_eIF3d_plusAux_plusHarr_4h_rep3_R1"
    "ribo_70_eIF4E_minusAux_plusHarr_4h_rep1_R1"
    "ribo_71_eIF4E_minusAux_plusHarr_4h_rep2_R1"
    "ribo_72_eIF4E_minusAux_plusHarr_4h_rep3_R1"
    "ribo_73_eIF4E_plusAux_plusHarr_4h_rep1_R1"
    "ribo_74_eIF4E_plusAux_plusHarr_4h_rep2_R1"
    "ribo_75_eIF4E_plusAux_plusHarr_4h_rep3_R1"
    "ribo_76_eIF4G1_plusAux_plusHarr_4h_rep1_R1"
    "ribo_77_eIF4G1_plusAux_plusHarr_4h_rep2_R1"
    "ribo_78_eIF4G1_plusAux_plusHarr_4h_rep3_R1"
    "ribo_79_eIF4G2_plusAux_plusHarr_4h_rep1_R1"
    "ribo_80_eIF4G2_plusAux_plusHarr_4h_rep2_R1"
    "ribo_81_eIF4G2_plusAux_plusHarr_4h_rep3_R1"
)

echo "Processing reads..."

# Loop through each sample to extract multimapped reads
for sample in "${samples[@]}"; do
    # Input and output file paths
    input_bam="results/split_bam/transcriptome/human/${sample}.bam"
    output_bam="results/multimapped_bam/transcriptome/human/multimapped_${sample}.bam"

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
    input_bam="results/split_bam/transcriptome/human/${sample}.bam"
    output_bam="results/unique_bam/transcriptome/human/unique_${sample}.bam"

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
    input_bam="results/split_bam/transcriptome/human/${sample}.bam"
    output_bam="results/primary_bam/transcriptome/human/primary_${sample}.bam"

    echo "Processing ${input_file} to extract primary reads..."

    # Extract primary alignments (flag 256 not set) and save to a new BAM file
    samtools view -h -b -F 256 "$input_bam" > "$output_bam"

    # Index the output BAM file
    samtools index "$output_bam"

    echo "Finished extracting primary reads from ${sample}. Output written to ${output_file}"
done

echo "Processing complete."
