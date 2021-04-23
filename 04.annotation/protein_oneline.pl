#!/usr/bin/perl
use strict;
my ($h, $s);
while (<STDIN>) {
	chomp;
	if (/^>/) {
		print STDOUT "$h\n$s\n" if $s;
		$h = $_;
		$s = "";
	} else {
		s/\*// if /\*$/;
		$s .= $_;
	}
}
print STDOUT "$h\n$s\n" if $s;
close STDIN;
