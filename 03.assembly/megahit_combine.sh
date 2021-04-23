#!/bin/bash
#$ --no-unal -S /bin/bash
cd /vol3/agis/fanwei_group/wangsen/projects/wheat_microbiome/rhizosphere_soil_metagenome/04.assembly/
cat Y*/*_R1_unmap.fq >> Combine_R1_unmap.fq
cat Y*/*_R2_unmap.fq >> Combine_R2_unmap.fq
rm Y*/*.fq
megahit -1 Combine_R1_unmap.fq -2 Combine_R2_unmap.fq --presets meta-large -t 20 -m 0.8 -o ./Unmap_combine
cd Unmap_combine && rm -rf intermediate_contigs
time bowtie2-build --threads 20 final.contigs.fa index
mv ../Combine*.fq .
echo 'Combine_R1_unmap.fq Combine_R2_unmap.fq'
time bowtie2 -p 20 -x index -1 Combine_R1_unmap.fq -2 Combine_R2_unmap.fq --local -k 1 --no-head --no-sq --no-unal -S /dev/null
rm index.* *.fq
