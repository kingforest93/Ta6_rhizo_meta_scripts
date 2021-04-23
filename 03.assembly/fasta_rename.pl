#!/usr/bin/perl
use strict;

defined $ARGV[0] || die "Usage: fasta_rename.pl input.fa\n";
open(IN, "<$ARGV[0]") || die "Cannot open file $ARGV[0]!\n";
my $num;
my %seqs;
while (<IN>) {
	chomp;
	if (/^>(\w+)/) {
		$num++;
		print ">Contig_$num\n";
	} else {
		print "$_\n";
	}
}
close(IN);
