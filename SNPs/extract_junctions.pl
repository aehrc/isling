#!/usr/bin/perl
# It's not a bug, it's a feature

use strict;
use warnings;

use Getopt::Long;

my $sam; #input sam file for alignment
my $ints; # file with integration sites
my $junctSam; #output sam file for juctions aligment

my $help;

GetOptions(	'ints=s' => \$ints,
		'sam=s' => \$sam,
		'junctSam=s'  => \$junctSam,
	        'help',    => \$help);

if ($help) { printHelp(); }

unless ($sam and $ints and $junctSam) { printHelp(); }

## first read integration site file
## extract the IDs of the reads and their sequences to grab from the alignment files
## store these in a hash in the format IDxxxseq

my %juncts;

open (INTS, $ints) || die "Could not open integration site file: $ints\n";
while (my $int = <INTS>) {

	chomp($int); #remove trailing newline

	#get fields from line
	my @parts = split("\t", $int);
	my $key = join("xxx", $parts[20], $parts[21]);
	

	#add to array of queries
	$juncts{$key} = 1;

}

close(INTS);

#check that there were some integrations

if ((scalar keys %juncts) == 0) { print "Warning: didn't find any integrations in $ints\n"; }


#check for matches between array of IDs and sequences, and an alignment file
#output new alignment file of matches

open (IN, $sam) || die "Could not open input host alignment file: $sam\n";
open (OUT, ">", $junctSam) || die "Could not open output host alingment file: $junctSam\n";
while (my $line = <IN>) {
	
	#output header lines
	if ($line =~ /^@/) { print OUT $line; next; }

	#get fields from line
	my ($flag, $ID, $seq) = (split("\t", $line))[1,0,9];

	if ($flag & 0x800) { next; } # skip supplementary alignments

	if ($flag & 0x4) { next; } # skip unmapped alignments
	
	#if sequence is reverse, reverse complement it
	if ($flag & 0x10) { ($seq) = reverseComp($seq); }
	
	#if IDxxxseq is in the hash, write the line
	my $key = join("xxx", $ID, $seq);
	if (exists $juncts{$key}) { print OUT $line; }

}
close(IN);
close(OUT);

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

sub reverseComp {
### return reverse compliment of sequence
	my ($seq) = @_;
	my $rev = reverse $seq;
	$rev =~ tr/ATGCatgc/TACGtacg/;
	return($rev);
}

