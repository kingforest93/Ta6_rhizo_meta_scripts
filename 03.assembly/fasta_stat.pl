#!/usr/bin/perl
use strict;

defined $ARGV[0] || die "Usage: fasta_stat.pl sequence.fasta\n";
open IN, "<", "$ARGV[0]" || die "Cannot open file $ARGV[0]!\n";
my ($seq_num, $seq_len);
my @seqs;
my $min = 10000; # give a default minimum length
my $max = 1; # give a default maximum length
while (<IN>) {
	chomp;
	if (/^>/) {
		$seq_num ++;
	} else {
		my $len = length $_;
		$seq_len += $len;
		push @seqs, $len;
		$min = $len if $len <= $min;
		$max = $len if $len >= $max;
	}
}
close IN;

print "Total number of sequences: $seq_num\n";
print "Total length of sequences: $seq_len bp\n";
print "Shortest sequence: $min bp\n";
print "Longest sequence: $max bp\n";
printf "Average sequence length: %d bp\n", $seq_len / $seq_num;

my ($Nn, $Ln, $cum);
@seqs = sort {$b <=> $a} @seqs;
foreach (@seqs) {
	$Ln += 1;
	$cum += $_;
	$Nn = $_;
	if ($cum >= 0.5 * $seq_len) {
		printf "N50 (sequence length): %d bp\n", $Nn;
		printf "L50 (the ith sequence): %d\n", $Ln;
		last;
	}
}
