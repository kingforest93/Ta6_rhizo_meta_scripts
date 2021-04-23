#!/usr/bin/perl
use strict;
defined $ARGV[0] || die "Usage: perl species2function_by_gene.pl gene.taxonomy gene.function\n";

print STDERR "parsing $ARGV[0] ...\n";
my %species;
open IN, "<", $ARGV[0] || die "Cannot open $ARGV[0]!\n";
while (<IN>) {
        chomp;
        next if /\|$/;
        my ($a, $b, $c) = split /\t/, $_;
        $c =~ s/.*\|//;
        push @{$species{$c}}, $a;
}
close IN;
print STDERR "done\n";

print STDERR "parsing $ARGV[1] ...\n";
my %function;
open IN, "<", $ARGV[1] || die "Cannot open $ARGV[1]!\n";
while (<IN>) {
	chomp;
	my ($a, $b) = split /\t/, $_;
	my @f = split /,/, $b;
	foreach my $f (@f) {
		$function{$a} = $f;
	}
}
close IN;
print STDERR "done\n";

print STDERR "link species to function ...\n";
foreach my $s (keys %species) {
	my %function_count;
	my $total_count = @{$species{$s}};
	my @function_list;
	foreach my $g (@{$species{$s}}) {
		$function_count{$function{$g}} += 1 if $function{$g};
	}
	foreach my $f (sort {$function_count{$b} <=> $function_count{$a}} keys %function_count) {
		push @function_list, "$function_count{$f}|$f";
	}
	print "$s\t$total_count\t" . join(",", @function_list) . "\n" if @function_list > 0;
}
print STDERR "done\n";
