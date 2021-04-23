#!/bin/bash
#$ -S /bin/bash
cd /vol3/agis/fanwei_group/wangsen/projects/wheat_microbiome/rhizosphere_soil_metagenome/05.annotation
DB=/public/agis/fanwei_group/wangsen/database/eggNOG
TMP=/dev/shm
#echo 'split the gene set into many parts, each with 1000000 sequences'
#time zcat final_gene_42_rep_annotations/sequence/final_gene_42_rep_part*.fasta.gz | split -a 2 -d -l 2000000 - part_42_
#for f in part_42_[0-9][0-9]; do
#	echo $f
#	diamond blastx -d ${DB}/eggnog_proteins.dmnd -q $f -p 72 -b 10 -c 1 -t ${TMP} --more-sensitive -e 0.001 --top 3 --query-cover 0 --subject-cover 0 -f 6 qseqid sseqid evalue bitscore -o $f.tbl
#	perl gene_best_hit.pl $f.tbl >> final_gene_42_rep.seed_orthologs
#	rm $f $f.tbl
#done
echo 'copy the eggNOG database to SSD'
cp ${DB}/eggnog.db ${TMP}/
emapper.py --data_dir ${TMP}/ --annotate_hits_table final_gene_42_rep.seed_orthologs --cpu 16 -o final_gene_42_rep
#rm ${TMP}/eggnog.db
