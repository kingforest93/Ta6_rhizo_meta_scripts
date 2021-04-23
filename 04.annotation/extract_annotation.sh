#!/bin/bash
cd /public/agis/fanwei_group/wangsen/projects/wheat_microbiome/rhizosphere_soil_metagenome/05.annotation
#TAXONKIT_DB=/public/agis/fanwei_group/wangsen/database/taxonomy
#echo "taxonomy lineage"
#time awk -F"\t" '$2!=""{gsub(/\..*/,"",$2);printf("%s\t%s\n",$1,$2)}' final_gene_42_rep.emapper.annotations | sed '1d' > taxid.list
#time cut -f2 taxid.list | taxonkit lineage | taxonkit reformat | cut -f3 | sed 's/;/|/g' > lineage.list
#time paste taxid.list lineage.list > final_gene_42_rep.taxonomy
#time perl -F"\t" -lane '$count{$F[2]}++; END{foreach(keys %count){$tax = $_; s/\|+/\t/g; print "$count{$tax}\t$_"}}' final_gene_42_rep.taxonomy > final_gene_42_rep.taxonomy.krona
#rm taxid.list lineage.list
#total=$(grep -c "^>" final_gene_42_rep.fasta)
#classfy=$(grep -c "^C" final_gene_42_rep.taxonomy)
#num=$(echo "$total - $classfy" | bc)
#sed -i "1i$num\tUnclassfied" final_gene_42_rep.taxonomy.krona
echo "KEGG KO"
time awk -F"\t" '$9!=""{printf("%s\t%s\n",$1,$9)}' final_gene_42_rep.emapper.annotations | sed '1d;s/ko://g' > final_gene_42_rep.KEGG_KO
echo "KEGG PATHWAY"
time awk -F"\t" '$10!=""{printf("%s\t%s\n",$1,$10)}' final_gene_42_rep.emapper.annotations | sed '1d;s/,map.*//g' > final_gene_42_rep.KEGG_PATHWAY
echo "KEGG BRITE"
time awk -F"\t" '$14!=""{printf("%s\t%s\n",$1,$14)}' final_gene_42_rep.emapper.annotations | sed '1d;s/br[0-9]\+,//g' > final_gene_42_rep.KEGG_BRITE
echo "eggNOG OG"
time awk -F"\t" '$19!=""{printf("%s\t%s\n",$1,$19)}' final_gene_42_rep.emapper.annotations | sed '1d;s/\@[0-9]\+//g' > temp.list
time perl -ne 'chomp; @F = split /[\t,]/, $_; $g = shift @F; @N = (); for $i (@F) {if ($i =~ /^\d/) {push @N, "ENOG50" . $i} else {push @N, $i}}; printf("%s\t%s\n", $g, join(",", @N))' temp.list > final_gene_42_rep.eggNOG_OG
rm temp.list
echo "COG GROUP"
time awk -F"\t" '$21!=""{printf("%s\t%s\t%s\n",$1,$21,$22)}' final_gene_42_rep.emapper.annotations | sed '1d' > temp.COG
time perl -F"\t" -lane '@c = split //, $F[1]; $c = join ",", @c; print "$F[0]\t$c\t$F[2]"' temp.COG > final_gene_42_rep.COG_GROUP
rm temp.COG
#echo "Length"
#time perl -ne 'if (/^>(\w+) #/) {$h=$1} else {$l = length($_); print "$h\t$l\n"}' final_gene_42_rep.fasta > final_gene_42_rep.length
echo "link species to function"
time perl species2function_by_gene.pl final_gene_42_rep.kaiju.taxonomy final_gene_42_rep.KEGG_KO > final_gene_42_rep.species2KEGG_KO
time perl species2function_by_gene.pl final_gene_42_rep.kaiju.taxonomy final_gene_42_rep.eggNOG_OG > final_gene_42_rep.species2eggNOG_OG
time perl species2function_by_gene.pl final_gene_42_rep.kaiju.taxonomy final_gene_42_rep.KEGG_BRITE > final_gene_42_rep.species2KEGG_BRITE
time perl species2function_by_gene.pl final_gene_42_rep.kaiju.taxonomy final_gene_42_rep.KEGG_PATHWAY > final_gene_42_rep.species2KEGG_PATHWAY
time perl species2function_by_gene.pl final_gene_42_rep.kaiju.taxonomy final_gene_42_rep.COG_GROUP > final_gene_42_rep.species2COG_GROUP
echo "link function to species"
time perl function2species_by_gene.pl final_gene_42_rep.KEGG_KO final_gene_42_rep.kaiju.taxonomy > final_gene_42_rep.KEGG_KO2species
time perl function2species_by_gene.pl final_gene_42_rep.eggNOG_OG final_gene_42_rep.kaiju.taxonomy > final_gene_42_rep.eggNOG_OG2species
time perl function2species_by_gene.pl final_gene_42_rep.KEGG_BRITE final_gene_42_rep.kaiju.taxonomy > final_gene_42_rep.KEGG_BRITE2species
time perl function2species_by_gene.pl final_gene_42_rep.KEGG_BRITE final_gene_42_rep.kaiju.taxonomy > final_gene_42_rep.KEGG_BRITE2species
time perl function2species_by_gene.pl final_gene_42_rep.COG_GROUP final_gene_42_rep.kaiju.taxonomy > final_gene_42_rep.COG_GROUP2species
