#!/usr/bin/perl
use strict;

die "Usage: perl abundance2taxonomy.pl gene2taxonomy.map sample1.abundance sample2.abundance ... \n" if not defined $ARGV[0];
my $f = shift @ARGV;
my %gene2tax;
print STDERR "parsing $f ...\n";
open IN, "<", $f or die "Cannot open $f!\n";
while (<IN>) {
	chomp;
	my ($a, $b, $c) = split /\t/, $_;
	@{$gene2tax{$a}} = split /\|/, $c;
}
close IN;
print STDERR "done\n";

foreach $f (@ARGV) {
	my %tax_abund;
	my $total;
	print STDERR "processing $f ... \n";
	open IN, "<", $f or die "Cannot open $f!";
	my @rank = ('d', 'p', 'c', 'o', 'f', 'g', 's');
	while (<IN>) {
		chomp;
		my ($a, $b) = split /\t/, $_;
		next if not $gene2tax{$a};
		next if $b eq "";
		my $r = 0;
		my $p = '';
		$total += $b;
		foreach my $t (@{$gene2tax{$a}}) {
			if ($t) {
				$p .= "$rank[$r]__$t|";
				$tax_abund{$p} += $b;
			}
			$r ++;
		}
	}
	close IN;
	print STDERR "done\n";
	
	print STDERR "output taxonomy abundances ... \n";
	$f =~ s/abundance/tax_abundance/;
	open OUT, ">", $f or die "Cannot write to $f!\n";
	foreach my $t (sort keys %tax_abund) {
		my $abund = $tax_abund{$t} * 100 / $total;
		$t =~ s/\|$//g;
		printf OUT "%s\t%.6f\n", $t, $abund;
	}
	close OUT;
	print STDERR "done\n";
}
