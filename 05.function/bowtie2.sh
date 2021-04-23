#!/bin/bash
#$ -S /bin/bash
cd /vol3/agis/fanwei_group/wangsen/projects/wheat_microbiome/rhizosphere_soil_metagenome/06.function
echo 'mapping reads to genes'
#time bowtie2-build --large-index --threads 12 final_gene_42_rep.fasta final_gene_42_rep
for i in Y*_R1_nohost.fq.gz; do
	if [ -s ${i%_R1_nohost.fq.gz}_42.sam ]; then
		continue;
	fi
	echo $i ${i/_R1/_R2}
	time bowtie2 -p 10 -x final_gene_42_rep -1 $i -2 ${i/_R1/_R2} --local --no-head --no-unal -S ${i%_R1_nohost.fq.gz}_42.sam
done
time perl sam2count.pl final_gene_42_rep.length Y*_42.sam
rm Y*_42.sam
