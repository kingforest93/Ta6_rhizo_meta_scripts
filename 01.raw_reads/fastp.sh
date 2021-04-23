#!/bin/bash
#$ -S /bin/bash
cd /vol3/agis/fanwei_group/wangsen/projects/wheat_microbiome/rhizosphere_soil_metagenome/01.raw_reads
for f in *R1.fastq.gz
do
	fastp -i $f -o ${f/L[0-9]*R1.fastq/R1_clean.fq} -I ${f/R1/R2} -O ${f/L[0-9]*R1.fastq/R2_clean.fq} -Q -5 --cut_front_window_size 1 -r -W 4 -M 30 -c --adapter_fasta NEB_adapters.fa -l 36 -w 6 -V -j ${f%_L[0-9]*.gz}.json -h ${f%_L[0-9]*.gz}.html
done
mkdir fastp_out
mv *.html *.json fastp_out
