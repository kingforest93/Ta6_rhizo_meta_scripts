#!/usr/bin/perl
use strict;
defined $ARGV[0] || die "Usage: perl combine_abundance.pl header sample1.abundance sample2.abundance sample3.abundance ... \n";
my %abundance;
my $header = shift @ARGV;
my @column;

foreach my $i (@ARGV) {
	print STDERR "parsing the abundance table of $i ... \r";
	open IN, "<", $i || die "Cannot open file $i!\n";
        my $h = $i;
        $h =~ s/_42.*//;
	push @column, $h;
	while (<IN>) {
		chomp;
		my ($a, $b) = split "\t", $_;
		$abundance{$a}{$h} = $b;
	}
	close IN;
}
print STDERR "\ndone\n";

print STDERR "combining the abundance of each file ... \n";
my $col = join "\t", @column;
print STDOUT "$header\t$col\n";
foreach my $k (sort keys %abundance) {
	my @v;
	foreach my $c (@column) {
		if (not defined $abundance{$k}{$c}) {
			push @v, 0.000000;
		}
		else {
			push @v, $abundance{$k}{$c};
		}
	}
	my $v = join "\t", @v;
	print STDOUT "$k\t$v\n";
}
print STDERR "done\n";
