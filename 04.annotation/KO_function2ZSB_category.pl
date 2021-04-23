#!/usr/bin/perl
use strict;
defined $ARGV[0] or die "Usage: perl KO_function2ZSB_category.pl genome.species2KEGG_KO\n";
my (%Zn_trans, %Zn_ligand, %Zn_uptake);
my @Zn_trans = split /\t/, "K09815	K09816	K09817	K19971	K19972	K10830	K14714	K14720	K14710	K14688	K14689	K14690	K14691	K14692	K14693	K14695	K14696	K14697	K14709	K14711	K14712	K14713	K14715	K14716	K14717	K14718	K14719	K16074	K07238	K09823	K02076";
my @Zn_ligand = split /\t/, "K08197	K08225	K16090	K14271	K22374	K05953	K01647	K01648	K15230	K15231	K05942	K01682	K00031	K00030	K17753	K01655	K02594	K10977	K01637	K01659	K01720	K20455	K03417	K10978	K05824	K24108	K23372	K23446	K23181	K23182	K23183	K23184	K01910	K01646	K01644	K01643	K09477	K11639	K15100	K03288	K16091	K00614	K00906	K11390	K13795	K23376	K03300	K00025	K00026	K00024	K00116	K00027	K00028	K00029	K00051	K01638	K01649	K08692	K14067	K07246	K14471	K14472	K09011	K01703	K01704	K00052	K01702	K21359	K21360	K17748	K17749	K17750	K21934	K12910	K20930	K20454	K15741	K15742	K11616	K24180	K00754	K18313	K01463	K22135	K01569	K22133	K08177	K19070	K19071	K19072	K18702";
my @Zn_uptake = split /\t/, "K13946	K14484	K14486	K14487	K13947	K24139	K20772	K01762	K05933	K01505	K02586	K02591	K02588	K00531	K22896	K22897	K22898	K22899	K02592	K02587	K22903	K02595	K15861	K04488	K15790	K02596	K02593	K02597	K02585	K23916";
foreach (@Zn_trans) {
	$Zn_trans{$_} = 1;
}
foreach (@Zn_ligand) {
	$Zn_ligand{$_} = 1;
}
foreach (@Zn_uptake) {
	$Zn_uptake{$_} = 1;
}
open IN, "<$ARGV[0]" or die "Cannot open $ARGV[0]!\n";
printf STDOUT "%-50s\t%20s\t%20s\t%20s\t%20s\n", "species_scientific_name","Total_genes","Zn_transport","Zn_ligand_metabolism","Zn_uptake_promotion";
while (<IN>) {
	chomp;
	my ($sp, $tot, $Zn_trans, $Zn_ligand, $Zn_uptake);
	my @col = split /[\t,]/, $_;
	$sp = shift @col;
	$tot = shift @col;
	foreach (@col) {
		my ($num, $ko) = split /\|/, $_;
		if ($Zn_trans{$ko}) {
			$Zn_trans += $num;
		} elsif ($Zn_ligand{$ko}) {
			$Zn_ligand += $num;
		} elsif ($Zn_uptake{$ko}) {
			$Zn_uptake += $num;
		} else {
			next;
		}
	}
	printf STDOUT "%-50s\t%20d\t%20d\t%20d\t%20d\n", $sp, $tot, $Zn_trans, $Zn_ligand, $Zn_uptake;
}
close IN;
