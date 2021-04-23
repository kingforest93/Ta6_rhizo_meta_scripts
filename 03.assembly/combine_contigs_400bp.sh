#!/bin/bash
#$ -S /bin/bash
cd /vol3/agis/fanwei_group/wangsen/projects/wheat_microbiome/rhizosphere_soil_metagenome/04.assembly/
echo 'mapping reads to combined contigs'
#time bowtie2-build --large-index --threads 20 combined_contigs_400bp.fa combined_contigs_400bp
#for i in Y*_R1_nohost.fq.gz; do
#	if [ -s ${i%_R1_nohost.fq.gz}.sam ]; then
#		continue
#	fi
#	echo $i ${i/_R1/_R2}
#	time bowtie2 -p 16 -x combined_contigs_400bp -1 $i -2 ${i/_R1/_R2} --local --no-unal --no-head -S | cut -f3,4,6 > ${i%_R1_nohost.fq.gz}.sam
#done
time perl sam2coverage.pl combined_contigs_400bp.length Y*.sam
#time perl combine_coverage.pl combined_contigs_400bp.length Y*.coverage > combined_400bp.coverage
paste Y*.coverage | cut -f 1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112 > temp.coverage
perl sum_coverage.pl combined_contigs_400bp.length temp.coverage > combined_400bp.coverage
rm temp.coverage
