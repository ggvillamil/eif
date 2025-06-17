#! /bin/bash

#SBATCH --job-name=mapping_riboseq_eif3d001
#SBATCH --mail-user=gabriel.villamil@mdc-berlin.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=logs/snakemake/%j.snakemake.stdout.log
#SBATCH --export=ALL
#SBATCH --chdir=.
#SBATCH --ntasks=72
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=100000MB


# Source .bashrc file
source ~/.bashrc

# Activate Snakemake conda environment
conda activate z_snakemake

# Run Snakemake
snakemake -j 72 -k -p --restart-times 1 --max-jobs-per-second 5 --rerun-incomplete --use-singularity

# Call another bash script
# conda activate ribopipe
# bash src/extract_multimappers_rnaseq.sh
