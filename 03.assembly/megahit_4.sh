#!/bin/bash
#$ -S /bin/bash

cd /vol3/agis/fanwei_group/wangsen/projects/wheat_microbiome/rhizosphere_soil_metagenome/04.assembly
for i in *[L]-4-*_R1_nohost.fq.gz; do
	echo $i ${i/R1_/R2_}
	megahit -1 $i -2 ${i/R1_/R2_} -o ${i%_R1_nohost.fq.gz} --presets meta-large -t 10 -m 0.1 --continue
	cd ${i%_R1_nohost.fq.gz} && rm -rf intermediate_contigs
        time bowtie2-build --threads 20 final.contigs.fa index
        time bowtie2 -p 20 -x index -1 ../$i -2 ../${i/R1_/R2_} --local -k 1 --no-head --no-sq --no-unal -S /dev/null --un-conc temp
        rm index.*
        mv temp.1 ${i%_R1_nohost.fq.gz}_R1_unmap.fq
        mv temp.2 ${i%_R1_nohost.fq.gz}_R2_unmap.fq
        cd ..
done
