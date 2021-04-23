#!/usr/bin/perl
use strict;
defined $ARGV[0] || die "Usage: perl function2species_by_gene.pl gene.function gene.taxonomy\n";

print STDERR "parsing $ARGV[0] ...\n";
my %function;
open IN, "<", $ARGV[0] || die "Cannot open $ARGV[0]!\n";
while (<IN>) {
	chomp;
	my ($a, $b) = split /\t/, $_;
	my @f = split /,/, $b;
	foreach my $f (@f) {
		push @{$function{$f}}, $a;
	}
}
close IN;
print STDERR "done\n";

print STDERR "parsing $ARGV[1] ...\n";
my %species;
open IN, "<", $ARGV[1] || die "Cannot open $ARGV[1]!\n";
while (<IN>) {
	chomp;
	next if /\|$/;
	my ($a, $b, $c) = split /\t/, $_;
	$c =~ s/.*\|//;
	$species{$a} = $c;
}
close IN;
print STDERR "done\n";

print STDERR "link function to species ...\n";
foreach my $f (keys %function) {
	my %species_count;
	my $total_count = @{$function{$f}};
	my @species_list;
	foreach my $g (@{$function{$f}}) {
		$species_count{$species{$g}} += 1 if $species{$g};
	}
	foreach my $s (sort {$species_count{$b} <=> $species_count{$a}} keys %species_count) {
		push @species_list, "$species_count{$s}|$s";
	}
	print "$f\t$total_count\t" . join(",", @species_list) . "\n" if @species_list > 0;
}
print STDERR "done\n";
