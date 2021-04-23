#!/usr/bin/perl
use strict;
defined $ARGV[0] or die "Usage: perl OG_function2ZSB_category.pl OGs_ZSB.txt genome.species2KEGG_KO\n";
my (%Zn_trans, %Zn_ligand, %Zn_uptake);
open IN, "<$ARGV[0]" or die "Cannot open $ARGV[0]!\n";
while (<IN>) {
	chomp;
	my @t = split /\t/, $_;
	if ($t[1] =~ /transmembrane/) {
		$Zn_trans{$t[0]} = 1;
	} elsif ($t[1] =~ /metabolism/) {
		$Zn_ligand{$t[0]} = 1;
	} else {
		$Zn_uptake{$t[0]} = 1;
	}
}
close IN;
open IN, "<$ARGV[1]" or die "Cannot open $ARGV[1]!\n";
printf STDOUT "%-50s\t%20s\t%20s\t%20s\n", "species_scientific_name","Zn_transport","Zn_ligand_metabolism","Zn_uptake_promotion";
while (<IN>) {
	chomp;
	my ($sp, $Zn_trans, $Zn_ligand, $Zn_uptake);
	my @col = split /[\t,]/, $_;
	$sp = shift @col;
	shift @col;
	foreach (@col) {
		my ($num, $og) = split /\|/, $_;
		if ($Zn_trans{$og}) {
			$Zn_trans += $num;
		} elsif ($Zn_ligand{$og}) {
			$Zn_ligand += $num;
		} elsif ($Zn_uptake{$og}) {
			$Zn_uptake += $num;
		} else {
			next;
		}
	}
	printf STDOUT "%-50s\t%20d\t%20d\t%20d\n", $sp, $Zn_trans, $Zn_ligand, $Zn_uptake;
}
close IN;
