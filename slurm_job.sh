#! /bin/bash

#SBATCH --job-name=extract_multimappers_rnaseq
#SBATCH --mail-user=gabriel.villamil@mdc-berlin.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=logs/snakemake/%j.snakemake.stdout.log
#SBATCH --export=ALL
#SBATCH --chdir=.
#SBATCH --ntasks=1
#SBATCH --time=168:00:00
#SBATCH --mem-per-cpu=256000MB


# Source .bashrc file
source ~/.bashrc

# Activate Snakemake conda environment
# conda activate z_snakemake

# Run Snakemake
# snakemake -j 1 -k -p --restart-times 1 --max-jobs-per-second 5 --rerun-incomplete --use-singularity

# Call another bash script
conda activate ribopipe
bash src/extract_multimappers_rnaseq.sh
