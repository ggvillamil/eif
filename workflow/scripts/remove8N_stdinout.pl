#!/usr/bin/perl

#Removes 4nt randomers from each side of sequencing reads in fastq files
#Appends randomers to read name
#Takes one argument, the file name
#Written by Emanuel Wyler, emanuel.wyler@alumni.ethz.ch

#Modified to use STDIN as input and STDOUT as output
#Modified by Gabriel Villamil, gabriel.villamil@mdc-berlin.de

use strict;

my $N1;
my $N2;

while (my $line1 = <STDIN>) {
	my $line2 = <STDIN>;
	my $line3 = <STDIN>;
	my $line4 = <STDIN>;
	chomp $line1;
	chomp $line2;
	chomp $line3;
	chomp $line4;
	$N1=substr($line2, 0, 4);
	$N2=substr($line2, length($line2)-4, 4);
	$line2 = substr($line2, 4, length($line2)-8);
	$line4 = substr($line4, 4, length($line4)-8);
	$line1 = $line1 . "_" . "$N1:$N2";
	print STDOUT "$line1\n$line2\n$line3\n$line4\n";
}

