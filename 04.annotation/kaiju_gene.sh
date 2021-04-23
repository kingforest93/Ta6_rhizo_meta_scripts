#!/bin/bash
#$ -S /bin/bash
cd /vol3/agis/fanwei_group/wangsen/projects/wheat_microbiome/rhizosphere_soil_metagenome/05.annotation
export KAIJU_DB=/public/agis/fanwei_group/wangsen/database/kaiju_db
kaiju -z 16 -a greedy -e 3 -s 65 -v -t $KAIJU_DB/nodes.dmp -f $KAIJU_DB/nr.fmi -o temp.kaiju -i ./final_gene_21_rep/final_gene_21_rep.fasta
time awk -F"\t" '{num[$1]++;tot++} END{print "T\t" tot;for(i in num){printf("%s\t%d\t%.2f\n", i, num[i], num[i]/tot*100)}}' temp.kaiju
time kaiju-addTaxonNames -v -t $KAIJU_DB/nodes.dmp -n $KAIJU_DB/names.dmp -r superkingdom,phylum,class,order,family,genus,species -i temp.kaiju | cut -f2,3,8 > temp.kaiju_names
time sed 's/NA.*//g;s/; /|/g;s/|$//g' temp.kaiju_names > temp.kaiju_taxonomy
time perl -F"\t" -lane 'if($F[2]=~/^\w/){$count{$F[2]}++}else{$un++}; END{print "$un\tUnclassfied"; foreach(sort keys %count){$tax = $_; s/\|/\t/g; print "$count{$tax}\t$_"}}' temp.kaiju_taxonomy > final_gene_21_rep.kaiju.taxonomy.krona
time awk -F"\t" '$3~/^\w/{print $0}' temp.kaiju_taxonomy > final_gene_21_rep.kaiju.taxonomy
rm temp.kaiju*
