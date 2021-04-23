#!/bin/bash
#$ -S /bin/bash
cd /vol3/agis/fanwei_group/wangsen/projects/wheat_microbiome/rhizosphere_soil_metagenome/05.annotation
time split -a 2 -l 2000000 -d combined_contigs_400bp.fa part_400bp_
time find . -name "part_400bp_*" | xargs -P 48 -I {} prodigal -p meta -a {}.faa -d {}.fna -f gbk -o {}.gbk -i {}
cat part_400bp_*.fna >> temp_gene.fna
cat part_400bp_*.faa >> temp_protein.faa
rm part_400bp_*
time perl fasta_filter_oneline.pl temp_gene.fna 100 > combined_400bp_gene_200bp.fna
rm temp_gene.fna
time perl fasta_rep_select.pl combined_400bp_gene_200bp.fna temp_protein.faa > combined_400bp_protein_200bp.faa
rm temp_protein.faa
time perl gene_complete_stat.pl combined_400bp_gene_200bp.fna > combined_400bp_gene_200bp.fna.stat
time perl fasta_length_stat.pl combined_400bp_gene_200bp.fna >> combined_400bp_gene_200bp.fna.stat
