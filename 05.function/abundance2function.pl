#!/usr/bin/perl
use strict;

die "Usage: perl abundance2taxonomy.pl gene2taxonomy.map sample1.abundance sample2.abundance ... \n" if not defined $ARGV[0];
my $g = shift @ARGV;
my %gene2function;
print STDERR "parsing $g ...\n";
open IN, "<", $g or die "Cannot open $g!\n";
while (<IN>) {
	chomp;
	my ($a, $b) = split /\t/, $_;
	@{$gene2function{$a}} = split /,/, $b;
}
close IN;
print STDERR "done\n";
$g =~ s/final_gene_42_rep\.//;
$g .= "_abundance";

foreach my $f (@ARGV) {
	my (%function_abund, %function_count);
	my $total;
	print STDERR "processing $f ... \n";
	open IN, "<", $f or die "Cannot open $f!\n";
	while (<IN>) {
		chomp;
		my ($a, $b) = split /\t/, $_;
		next if not $gene2function{$a};
		next if $b eq "";
		foreach my $i (@{$gene2function{$a}}) {
			$function_abund{$i} += $b;
			$function_count{$i} += 1;
		}
	}
	close IN;
	print STDERR "done\n";
	
	foreach my $i (keys %function_abund) {
		$function_abund{$i} = $function_abund{$i} / $function_count{$i};
		$total += $function_abund{$i};
	}
	
	print STDERR "output function abundances ... \n";
	$f =~ s/abundance/$g/;
	open OUT, ">", $f or die "Cannot write to $f!\n";
	foreach my $i (keys %function_abund) {
		my $abund = $function_abund{$i} * 100 / $total;
		printf OUT "%s\t%.6f\n", $i, $abund;
	}
	close OUT;
	print STDERR "done\n";
}
