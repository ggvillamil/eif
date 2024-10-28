#! /bin/bash

#SBATCH --job-name=hct116_anno_orfquant_reports
#SBATCH --mail-user=gabriel.villamil@mdc-berlin.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=logs/snakemake/%j.snakemake.stdout.log
#SBATCH --export=ALL
#SBATCH --chdir=.
#SBATCH --ntasks=28
#SBATCH --time=168:00:00
#SBATCH --mem-per-cpu=10000MB


# Source .bashrc file
source ~/.bashrc

# Activate Snakemake conda environment
conda activate z_snakemake

# Run Snakemake
snakemake -j 28 -k -p --restart-times 1 --max-jobs-per-second 5 --rerun-incomplete --use-singularity
