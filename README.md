# eIF
Analysis pipeline for the processing and analysis of ribosome profiling and total RNA sequencing data from degron cells targeting eukaryotic translation initiation factors (eIFs). It is implemented on Snakemake and is executed in a [Docker](https://www.docker.com/resources/what-container/) container for distribution and reproducibility.

## Features
The workflow takes in ribosome profiling reads in FASTQ format as input and automates read processing, mapping, and quality-checking. Mapped reads are analyzed with RiboseQC, ORFquant, and Ribostan to calculate P-site offsets, predict open reading frames (ORFs) on transcripts, including upstream ORFs (uORFs), quantify ribosome occupancy, and with total RNA-seq data, calculate differential translation efficiency.


### Requirements
1. [Snakemake](https://snakemake.readthedocs.io/en/stable/)
2. [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) or Docker


# Install
### Clone the Repo
This project is structured as a self-contained repository. Start by cloning the repository with this command:

	git clone git@github.com:ggvillamil/eif.git

### Download raw data

# Use
### Run the Workflow
#### Submit job

	sbatch slurm_job.sh

#### Run in an interactive session

	srun --pty bash
	cd eif
	conda activate snakemake
	snakemake --use-singularity --cores 16

### Run by Rule

	srun --pty bash
	cd eif
	conda activate snakemake
	snakemake --use-singularity --cores 16

# Contributing

