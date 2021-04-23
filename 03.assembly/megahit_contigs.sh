#!/bin/bash
#$ -S /bin/bash
cd /vol3/agis/fanwei_group/wangsen/projects/wheat_microbiome/rhizosphere_soil_metagenome/04.assembly/
megahit -r combined_contigs_200bp.fa --presets meta-large -t 36 -m 0.3 -o ./assembly_contigs
#mv assembly_contigs/final.contigs.fa assembly_contigs_200bp.fa
#rm -rf combine_total
time perl fasta_stat.pl assembly_contigs_200bp.fa > assembly_contigs_200bp.fa.stat
time perl fasta_filter.pl assembly_contigs_200bp.fa > assembly_contigs_400bp.fa
time perl fasta_stat.pl assembly_contigs_400bp.fa > assembly_contigs_400bp.fa.stat
#echo 'mapping reads to combined contigs'
#time bowtie2-build --large-index --threads 12 assembly_contigs_200bp.fa index
#for i in YS-2*R1_nohost.fq.gz; do
#	echo $i ${i/_1/_2}
#	time bowtie2 -p 24 -x index -1 $i -2 ${i/_R1/_R2} --local --very-fast-local -k 1 --no-head --no-sq --no-unal -S /dev/null
	#time bowtie2 -p 24 -x index -1 $i -2 ${i/_R1/_R2} --local --no-unal | samtools sort -@ 24 -o ${i%_R1_nohost.fq.gz}.bam
#done
#rm index.*
