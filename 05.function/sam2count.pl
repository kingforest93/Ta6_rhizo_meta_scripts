#!/usr/bin/perl
use strict;

die "Usage: perl sam2count.pl gene2length.tsv sample1.sam sample2.sam ... \n" if not defined $ARGV[0];

my $f = shift @ARGV;
print STDERR "processing $f ... \n";
open IN, "<", $f or die "Cannot open $f!\n";
my %gene_length;
while (<IN>) {
	chomp;
	my ($a, $b) = split /\t/, $_;
	$gene_length{$a} = $b;
}
close IN;
print STDERR "done\n";

foreach $f (@ARGV) {
	my %gene_count;
	my %single_unmap = ('73' => '1', '133' => '1', '89' => '1', '121' => '1', '165' => '1', '181' => '1', '101' => '1', '117' => '1', '153' => '1', '185' => '1', '69' => '1', '137' => '1');
	#my %double_unmap = ('77' => '1', '141' => '1');
	print STDERR "processing $f ... \n";
	open IN, "<", $f or die "Cannot open $f!";
	while (<IN>) {
		chomp;
		my @F = split /\t/, $_;
		#next if exists $double_unmap{$F[1]};
		#my ($mis, $match, $indel, $align, $identity);
		#$mis = ($F[11] =~ /XM:i:(\d+)/);
		#while ($F[5] =~ /(\d+)M/g) {$match += $1};
		#$match = $match - $mis;
		#while ($F[5] =~ /(\d+)[DI]/g) {$indel += $1};
		#$align = $match + $mis + $indel;
		#$identity = $match / $align;
		#next if ($align <= 45 || $identity <= 0.95);
		if (exists $single_unmap{$F[1]}) {
			$gene_count{$F[2]} += 1;
		} elsif ($F[6] ne "=") {
			$gene_count{$F[2]} += 1;
		} else {
			$gene_count{$F[2]} += 0.5;
		}
	}
	close IN;
	print STDERR "done\n";
	
	print STDERR "output gene counts ... \n";
	$f =~ s/sam/count/;
	open OUT, ">", $f or die "Cannot open $f!";
	foreach my $g (keys %gene_count) {
		print OUT "$g\t$gene_length{$g}\t$gene_count{$g}\n";
	}
	close OUT;
	print STDERR "done\n";
}
