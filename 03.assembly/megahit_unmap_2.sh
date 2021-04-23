#!/bin/bash
#$ -S /bin/bash
cd /vol3/agis/fanwei_group/wangsen/projects/wheat_microbiome/rhizosphere_soil_metagenome/04.assembly/
cat YL-2-*/*_R1_unmap.fq >> YL-2_R1.fq
cat YL-2-*/*_R2_unmap.fq >> YL-2_R2.fq
rm YL-2-*/*.fq
megahit -1 YL-2_R1.fq -2 YL-2_R2.fq --presets meta-large -t 20 -m 0.8 -o ./YL-2_unmap
cd YL-2_unmap && rm -rf intermediate_contigs
time bowtie2-build --threads 20 final.contigs.fa index
mv YL-2*.fq .
echo 'YL-2_R1.fq YL-2_R2.fq'
time bowtie2 -p 20 -x index -1 YL-2_R1.fq -2 YL-2_R2.fq --local -k 1 --no-head --no-sq --no-unal -S /dev/null --un-conc temp
mv temp.1 YL-2_R1_unmap.fq
mv temp.2 YL-2_R2_unmap.fq
rm index.*
