# Ribopipe2
Ribopipe2 is an automated workflow for processing and analyzing ribosome profiling data. It is implemented on Snakemake and is executed in a dedicated [Docker](https://www.docker.com/resources/what-container/) container for distribution and reproducibility.

## Features
The workflow takes in ribosome profiling reads in FASTQ format as input and automates read processing, mapping, and quality-checking. Mapped reads are further analyzed with RiboseQC, ORFquant, and RiboStan to calculate P-site offsets, predict open reading frames (ORFs) on transcripts, and quantify ribosome occupancy. 

Ribopipe2 also calculates translation efficiency of ORFs when bulk RNA-seq data is provided.


### Requirements
1. [Snakemake](https://snakemake.readthedocs.io/en/stable/) v7.28.3
2. [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html)


# Install
### Install Miniconda
Download and install [Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) with the following:

	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda_installer.sh
	bash miniconda_installer.sh -b


### Install Snakemake
The recommended way to install Snakemake is via Mamba. If your system already uses Conda but not Mamba run this first:

	conda install -n base -c conda-forge mamba

Then you can install Snakemake with:

	conda activate base
	mamba create -c conda-forge -c bioconda -n snakemake snakemake=7.28.3

More detailed information on Snakemake installation can be found [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).


### Clone the Repo
Ribopipe2 is structured as a self-contained repository for an independent ribosome profiling project, which may include multiple samples and replicates. The parent directory is designated as the project folder. Start by cloning the repository with this command:

	git clone https://github.com/ohlerlab/Ribopipe2.git <project-folder>


# Use
### Prepare the Sample Sheet
### Configure the Workflow
### Run the Workflow

#### Submit job

	qsub snake_job.sh

#### Run in an interactive job

	qrsh
	cd <path-to-project-folder>
	conda activate snakemake
	snakemake --use-singularity --cores 16

### Run by Rule

## Contributing


---


## Additional Notes

<details>
	<summary>Show</summary>

### Dockerfile
The Dockerfile used to build the Docker image is included in the `/docker` folder. You can edit the Dockerfile and re-build the Docker image to customize the container Ribopipe2 runs on.


### Step Into Container
You can step into Ribopipe2's container without executing the Snakemake worfklow with the following command.

	singularity run -B <path-to-project-folder>:<path-to-project-folder> docker://ggvillamil/ribopipe2

The `-B` flag binds your project folder so that files and folders within the project folder are accessible inside the container.

</details>

