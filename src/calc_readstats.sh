#!/bin/bash




rm -rf results/post/readstats.csv

echo "sample,raw_reads,cutadapt_reads,collapse_reads,trim_reads,filter_reads,reads_mapped" > results/post/readstats.csv


for sample in $(tail -n +2 config/samples.csv)
do
	var0=$(echo $(zcat data/riboseq/$sample.fastq.gz | wc -l)/4|bc)
	var1=$(echo $(zcat results/cutadapt_reads/$sample.fastq.gz | wc -l)/4|bc)
	var2=$(echo $(zcat results/collapse_reads/$sample.fastq.gz | wc -l)/4|bc)
	var3=$(echo $(zcat results/trim_reads/$sample.fastq.gz | wc -l)/4|bc)
	var4=$(echo $(zcat results/filter_reads/$sample/$sample.fastq.gz | wc -l)/4|bc)
	var5=$(awk $'/SN\treads mapped:/{ print $4 }' results/star/genome/reports/$sample/$sample.bamstats.txt)

	echo "$sample,$var0,$var1,$var2,$var3,$var4,$var5" >> results/post/readstats.csv
done
