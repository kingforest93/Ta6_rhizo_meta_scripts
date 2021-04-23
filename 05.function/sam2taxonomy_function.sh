#!/bin/bash
#$ -S /bin/bash
cd /vol3/agis/fanwei_group/wangsen/projects/wheat_microbiome/rhizosphere_soil_metagenome/06.function 
#time perl sam2count.pl final_gene_42_rep.length Y*_42.sam
#time perl combine_count.pl Gene Y*_42.count > combined_42.count
#time perl count2taxonomy.pl final_gene_42_rep.kaiju.taxonomy Y*_42.count
#time perl combine_count.pl Taxonomy Y*_42.tax_count > combined_42.tax_count
#time perl count2function.pl final_gene_42_rep.KEGG_KO Y*_42.count
#time perl combine_count.pl KEGG_KO Y*_42.KEGG_KO_count > combined_42.KEGG_KO_count
#time perl count2function.pl final_gene_42_rep.KEGG_PATHWAY Y*_42.count
#time perl combine_count.pl KEGG_PATHWAY Y*_42.KEGG_PATHWAY_count > combined_42.KEGG_PATHWAY_count
#time perl count2function.pl final_gene_42_rep.KEGG_BRITE Y*_42.count
#time perl combine_count.pl KEGG_BRITE Y*_42.KEGG_BRITE_count > combined_42.KEGG_BRITE_count
#time perl count2function.pl final_gene_42_rep.eggNOG_OG Y*_42.count
#time perl combine_count.pl eggNOG_OG Y*_42.eggNOG_OG_count > combined_42.eggNOG_OG_count
time perl count2function.pl final_gene_42_rep.COG_GROUP Y*_42.count
time perl combine_count.pl COG_GROUP Y*_42.COG_GROUP_count > combined_42.COG_GROUP_count
#for i in Y*_42.count; do
#	echo $i
#	time perl count2abundance.pl $i > ${i/count/abundance}
#done
#time perl combine_abundance.pl Gene Y*_42.abundance > combined_42.abundance
#time perl abundance2taxonomy.pl final_gene_42_rep.kaiju.taxonomy Y*_42.abundance
#time perl combine_abundance.pl Taxonomy Y*_42.tax_abundance > combined_42.tax_abundance
#time perl abundance2function.pl final_gene_42_rep.KEGG_KO Y*_42.abundance
#time perl combine_abundance.pl KEGG_KO Y*_42.KEGG_KO_abundance > combined_42.KEGG_KO_abundance
time perl abundance2function.pl final_gene_42_rep.KEGG_PATHWAY Y*_42.abundance
time perl combine_abundance.pl KEGG_PATHWAY Y*_42.KEGG_PATHWAY_abundance > combined_42.KEGG_PATHWAY_abundance
time perl abundance2function.pl final_gene_42_rep.KEGG_BRITE Y*_42.abundance
time perl combine_abundance.pl KEGG_BRITE Y*_42.KEGG_BRITE_abundance > combined_42.KEGG_BRITE_abundance
#time perl abundance2function.pl final_gene_42_rep.eggNOG_OG Y*_42.abundance
#time perl combine_abundance.pl eggNOG_OG Y*_42.eggNOG_OG_abundance > combined_42.eggNOG_OG_abundance
time perl abundance2function.pl final_gene_42_rep.COG_GROUP Y*_42.abundance
time perl combine_abundance.pl COG_GROUP Y*_42.COG_GROUP_abundance > combined_42.COG_GROUP_abundance
