#!/usr/bin/perl -w
use strict;
defined $ARGV[0] or die "Usage: perl read_stat.pl fastp.out.log\n";
my (@read1, @read2);
open IN, "<", $ARGV[0] or die "Cannot open $ARGV[0]!\n";
my ($r1b, $r2b, $r1a, $r2a);
print "Sample\tRead\tNum_reads_before\tNum_bases_before\tNum_reads_after\tNum_bases_after\n";
while (<IN>) {
	chomp;
	$r1b = 1 if /^Read1 before/;
	push @read1, $1 if /^total reads: (\d+)/ && $r1b;
	push @read1, $1 if /^total bases: (\d+)/ && $r1b;
	if (/^Read2 before/) {$r1b = 0; $r2b = 1};
        push @read2, $1 if /^total reads: (\d+)/ && $r2b;
        push @read2, $1 if /^total bases: (\d+)/ && $r2b;
	if (/^Read1 after/){$r2b = 0; $r1a = 1};
	push @read1, $1 if /^total reads: (\d+)/ && $r1a;
	push @read1, $1 if /^total bases: (\d+)/ && $r1a;
	if (/^Read2 after/) {$r1a = 0; $r2a = 1};	
	push @read2, $1 if /^total reads: (\d+)/ && $r2a;
	push @read2, $1 if /^total bases: (\d+)/ && $r2a;
	if (/^JSON report: (\S+)\.json/) {
		$r2a = 0;
		print "$1\tRead1\t" . join("\t", @read1) . "\n";
		print "$1\tRead2\t" . join("\t", @read2) . "\n";
		@read1 = @read2 = ();
	}
}
