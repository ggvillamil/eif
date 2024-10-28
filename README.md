# eIF
Analysis pipeline for the processing and analysis of ribosome profiling and total RNA sequencing data from degron cells targeting eukaryotic translation initiation factors (eIFs). It is implemented on Snakemake and is executed in a [Docker](https://www.docker.com/resources/what-container/) container for distribution and reproducibility.

### Features
The workflow takes in ribosome profiling reads in FASTQ format as input and automates read processing, mapping, and quality-checking. Mapped reads are analyzed with RiboseQC, ORFquant, and Ribostan to calculate P-site offsets, predict open reading frames (ORFs) on transcripts, including upstream ORFs (uORFs), quantify ribosome occupancy, and with total RNA-seq data, calculate differential translation efficiency.


### Requirements
1. [Snakemake](https://snakemake.readthedocs.io/en/stable/)
2. [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) or [Docker](https://docs.docker.com/get-started/docker-overview/)


## Install
### Clone the Repo
This project is structured as a self-contained repository. Start by cloning the repository with this command:

	git clone git@github.com:ggvillamil/eif.git

### Download raw data
Download raw ribosome profiling and RNA sequencing fastq files:

	cd eif
	wget http://bimsbstatic.mdc-berlin.de/ohler/gabriel/eif/data.tar.gz
	tar -zxvf data.tar.gz

### Download resource files
Download genome sequence, transcriptome sequence, transcript annotation, and other required resource files:

	cd eif
	wget http://bimsbstatic.mdc-berlin.de/ohler/gabriel/eif/resources.tar.gz
	tar -zxvf resources.tar.gz

## Use
### Run the Workflow
#### Submit job

	sbatch slurm_job.sh

#### Run in an interactive session

	srun --pty bash
	cd eif
	conda activate snakemake
	snakemake --use-singularity --cores 16

## To Do

- [ ] Use ribo-seq data to call translating ORFs and see global preference for uORFs/mORFs and start codon type
- [ ] Use harringtonine data to refine resolution and “validate” initiation sites (peak calling with ORF-RATER or re-implementation of ORF-RATER)
- [ ] Build predictive model to call TIS from harringtonine data
- [ ] Improve the current random forest model’s predictiveness (i.e. investigate two predicted populations)
- [ ] Include eIF4G2 and eIF4G3 in random forest models
- [x] Call uORFs from ribo-seq data (start immediately by using ORFquant on Slurm)
- [ ] Build model for 5’UTR-specific features
- [x] Pool ribo-seq libraries and run them through ORFquant
- [x] Use HCT116 transcript annotation from long read sequencing (if it looks good, adapt to the rest of the analyses)
- [ ] Compare responsive transcripts between 4- and 8-hour data sets in all factors
- [ ] Compare uORF/TIS changes  between 4- and 8-hour data sets in all factors