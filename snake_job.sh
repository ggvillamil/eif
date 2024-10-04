#! /bin/bash

#$ -N eif4f
#$ -M gabriel.villamil@mdc-berlin.de
#$ -m beas # send email beginning, end, and suspension
#$ -cwd #start from current directory
#$ -V #export all the environmental variables into the context of the job
#$ -j yes #merge the stderr with the stdout
#$ -o logs/snakemake/snakemake.stdout.log #stdout, job log
#$ -l h_vmem=20G
#$ -l h_rt=72:00:00
#$ -l os=centos7 -pe smp 2


source ~/.bashrc
conda activate snakemake
snakemake -j 64 -k -p --restart-times 1 --max-jobs-per-second 5 --rerun-incomplete --use-singularity --singularity-args "-B /fast/AG_Ohler/gabriel/eif4f/:/fast/AG_Ohler/gabriel/eif4f/" --cluster-config config/cluster.json --cluster="qsub -cwd -V -l m_mem_free={cluster.m_mem_free} -l h_rt={cluster.h_rt} -pe {cluster.pe} -j yes -o sge_log"