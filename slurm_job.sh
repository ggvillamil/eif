#! /bin/bash

#SBATCH --job-name=quantify_primary_riboseq_eif3d_markusfilter
#SBATCH --mail-user=gabriel.villamil@mdc-berlin.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=logs/snakemake/%j.snakemake.stdout.log
#SBATCH --export=ALL
#SBATCH --chdir=.
#SBATCH --ntasks=48
#SBATCH --time=168:00:00
#SBATCH --mem-per-cpu=30000MB


# Source .bashrc file
source ~/.bashrc

# Activate Snakemake conda environment
conda activate z_snakemake

# Run Snakemake
snakemake -j 48 -k -p --restart-times 1 --max-jobs-per-second 5 --rerun-incomplete --use-singularity

# Call another bash script
# conda activate ribopipe
# bash src/extract_multimappers_rnaseq.sh
