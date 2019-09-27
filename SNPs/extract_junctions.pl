#!/usr/bin/perl
# It's not a bug, it's a feature

use strict;
use warnings;

use Getopt::Long;

my $sam; #input sam file for alignment
my $ints; # file with integration sites
my $junctSam; #output sam file for juctions aligment

my $help;

GetOptions('ints=s' => \$ints,
		   'sam=s' => \$sam,
		   'junctSam=s'  => \$junctSam,
	           'help',    => \$help);

if ($help) { printHelp(); }

unless ($sam and $ints and $junctSam) { printHelp(); }

## first read integration site file
## extract the IDs of the reads and their sequences to grab from the alingment files
## store these in an array in the format IDxxxseq

my @juncts;

open (INTS, $ints) || die "Could not open integration site file: $ints\n";
while (my $int = <INTS>) {

	#get fields from line
	my @parts = split("\t", $int);

	#add to array of queries
	push(@juncts, join("xxx", $parts[18], $parts[19]));

}

close(INTS);

#then look through each aligment file (host and vector).  If there is a match between an element of @juncts and the line, output the line
checkAlignment($sam, $junctSam, @juncts);

exit;

#----------------------------------------------------------------------------------------------------------------------------------------
### Subroutines ###
#----------------------------------------------------------------------------------------------------------------------------------------

sub printHelp {
	print "Script for extraction of junction reads from host and viral alignments\n\n";
	print "Usage:\n";
	print "\tperl extract_junctions.pl --sam <sam> --ints <txt> --junctSam <out> --help\n\n";
	print "Arguments:\n";
	print "\t--sam:   Alignment to filter (sam)\n";
	print "\t--ints:  Integration site txt file\n";
	print "\t--junctSam:  Alingment containig only junction fragments (sam)\n";
	print "\t--help:    Print this help message and quit\n";

	exit;
}

sub checkAlignment{

#check for matches between array of IDs and sequences, and an alignment file
#output new alignment file of matches

my $infile = shift @_; #path to infile is first input
my $outfile = shift @_; #path to outfile is second input
my @queries = @_; #everything else is queries

open (IN, $infile) || die "Could not open input host alignment file: $infile\n";
open (OUT, ">", $outfile) || die "Could not open output host alingment file: $outfile\n";
while (my $line = <IN>) {
	
	#output header lines
	if ($line =~ /^@/) { print OUT $line; next; }

	#get fields from line
	my ($flag, $ID, $seq) = (split("\t", $line))[1,0,9];

	if ($flag & 0x800) { next; } # skip supplementary alignments

	#if read ID is a match
	my @matched = grep {/$ID/} @queries;
	
	if (@matched) {
		#if sequence is reverse, reverse complement it
		if ($flag & 0x10) { ($seq) = reverseComp($seq); }
		#if sequence is a match, write line
		if ( grep {/$seq/} @matched) { print OUT $line; }
	}
}
close(IN);
close(OUT);

}

sub reverseComp {
### return reverse compliment of sequence
	my ($seq) = @_;
	my $rev = reverse $seq;
	$rev =~ tr/ATGCatgc/TACGtacg/;
	return($rev);
}

