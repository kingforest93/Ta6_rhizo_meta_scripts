#!/bin/bash
##$ -S /bin/bash
cd /vol3/agis/fanwei_group/wangsen/projects/wheat_microbiome/rhizosphere_soil_metagenome/02.clean_reads
for i in *_R1_clean.fq.gz; do
	echo $i ${i/R1_/R2_}
	time bowtie2 -p 8 -x bowtie2_index/wheat_whole -1 $i -2 ${i/R1_/R2_} --local -k 1 --no-head --no-sq --no-mixed --no-unal -S /dev/null --un-conc-gz temp
	mv temp.1 ${i/clean/nohost}
	mv temp.2 ${i/R1_clean/R2_nohost}
	fq_stat.sh ${i/clean/nohost}
    fq_stat.sh ${i/R1_clean/R2_nohost}
done
