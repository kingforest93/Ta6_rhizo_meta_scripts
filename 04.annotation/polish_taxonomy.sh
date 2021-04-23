#!/bin/bash
cd /vol3/agis/fanwei_group/wangsen/projects/wheat_microbiome/rhizosphere_soil_metagenome/05.annotation/final_gene_42_rep_annotations/taxonomy
perl -lpne 's/\|{2,}.*//g; s/\|\w+\ssp\..*//g; s/\|(Candidatus|miscellaneous|uncultured|candidate)\s.*//g; s/\|\S+\s\S+\s\S+.*//g; s/\|$//g' final_gene_42_rep.kaiju.taxonomy > temp.taxonomy
time perl -F"\t" -lane 'if($F[2]=~/^\w/){$count{$F[2]}++}; END{foreach(sort keys %count){$tax = $_; s/\|/\t/g; print "$count{$tax}\t$_"}}' temp.taxonomy > temp.taxonomy.krona
mv temp.taxonomy final_gene_42_rep.kaiju.taxonomy
mv temp.taxonomy.krona final_gene_42_rep.kaiju.taxonomy.krona
