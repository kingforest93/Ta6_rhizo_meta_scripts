#!/bin/bash
#$ -S /bin/bash
cd /vol3/agis/fanwei_group/wangsen/projects/wheat_microbiome/rhizosphere_soil_metagenome/05.annotation
WORK=/vol3/agis/fanwei_group/wangsen/projects/wheat_microbiome/rhizosphere_soil_metagenome/05.annotation
TMP=/dev/shm
mmseqs createdb combined_400bp_gene_200bp.fna tmp_db --dbtype 2 --createdb-mode 0 --shuffle 0
mmseqs linclust tmp_db tmp_clu ${TMP}/temp --min-seq-id 0.95 -c 0.9 --cov-mode 1 --cluster-mode 2 --threads 24 --split-memory-limit 0 --force-reuse 0
mmseqs createsubdb tmp_clu tmp_db tmp_clu_rep
mmseqs convert2fasta tmp_clu_rep final_gene_42_rep.fasta
rm -rf ${TMP}/temp
time perl gene_complete_stat.pl final_gene_42_rep.fasta > final_gene_42_rep.fasta.stat
time perl fasta_rep_select.pl final_gene_42_rep.fasta combined_400bp_protein_200bp.faa > final_protein_42_rep.fasta
time perl fasta_length_stat.pl final_gene_42_rep.fasta >> final_gene_42_rep.fasta.stat
