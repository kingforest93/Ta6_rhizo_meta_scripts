#!/usr/bin/perl
use strict;

die "Usage: perl count2abund.pl sample.count\n" if not defined $ARGV[0];
my %gene_abund;
my $total_abund;

print STDERR "processing $ARGV[0] ... \n";
open IN, "<", "$ARGV[0]" or die "Cannot open $ARGV[0]!\n";
while (<IN>) {
	chomp;
	my ($a, $b, $c) = split /\t/, $_;
	$gene_abund{$a} = $c / $b;
	$total_abund += $c / $b;
}
close IN;
print STDERR "done\n";

print STDERR "finish and output gene abundances ... \n";
foreach (keys %gene_abund) {
	my $abund = $gene_abund{$_} * 100 / $total_abund;
	printf STDOUT "%s\t%.6f\n", $_, $abund;
}
print STDERR "done\n";
