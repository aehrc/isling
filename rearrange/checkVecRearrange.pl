#!/usr/bin/perl
# It's not a bug, it's a feature

use strict;
use warnings;

use ViralIntegration;
use Getopt::Long;

my $viral;
my $output = "rearrange.txt";
my $bed;
my $verbose;
my $thresh = 0.95;
my $help;

GetOptions(	   'thesh=f'  => \$thresh,
		   'viral=s'  => \$viral,
		   'output=s' => \$output,
		   'bed=s'	  => \$bed,
		   'verbose'  => \$verbose,
		   'help',    => \$help);

if ($help) { printHelp(); }

unless ($viral) { printHelp(); }

if (($thresh < 0) or ($thresh > 1)) { printHelp(); }


#open input file
open (VIRAL, $viral) || die "Could not open viral alignment file: $viral\n";
if ($verbose) { print "Processing viral alignment...\n"; }
my ($rearrange, $line, $bedline);

#open output file(s)
open (OUTFILE, ">$output") || die "Could not open output file: $output\n";
if ($bed) { open (BEDFILE, ">$bed") || die "Could not open output file: $output\n"; }

print OUTFILE "ReadID\tpossibleVecRearrange\tTotalGapBP\tGaps\tsups+secs\tseq\n";

#read input file and write to output file
while (my $vl = <VIRAL>) {
	if ($vl =~ /^@/) { next; } # skip header lines

	my @parts = split("\t", $vl);
	

	if ($parts[1] & 0x800) { next(); } # skip supplementary alignments

	if ($parts[2] eq "*") { next; } # skip unaligned reads

	my $cig = $parts[5]; 

	#get information about location of alignment
	my $ref = $parts[2];
	my $pos = $parts[3];

	#get read ID
	my $ID = $parts[0];
	
	#get supplementary (SA) and secondary (XA) alignments in order to check for possible vector rearrangements
	my ($sec, $sup) = getSecSup($vl);
	
	#get sequence and orientation
	my $seq;
	my $ori;
	#reverse complement sequence if necessary
	if ($parts[1] & 0x10) { 
		($seq) = reverseComp($parts[9]); 
		$ori   = 'r';
	}
	else				  { 
		$seq = $parts[9]; 
		$ori = 'f';
	}
	
	#check for possible rearrangements
	my @reData = isRearrange($cig, $ori, $ref, $pos, (join(";", $sup, $sec)), $seq, $thresh);
	my $isRearrange = shift @reData;
	my $gapBP = shift @reData;
	my $gaps = shift @reData;
	my @aligns = @reData;
		
	# print to bed file
	if ($bed) {			
		foreach my $align (@aligns) {
			#get info from this alignment
			my ($start, $stop, $ref, $pos, $sense, $cig) = split("xxx", $align);
				
			#get coordinates in virus
			my ($gStart, $gStop) = getGenomicCoords($start, $stop, $pos, $sense, $cig);
				
			#make line and print to file
			$bedline = join("\t", $ref, $gStart, $gStop, $isRearrange, $ID);
			print BEDFILE "$bedline\n";
		}		
	}

	#print to txt file
	$line = join("\t", $ID, $isRearrange, $gapBP, $gaps, $sup, $seq);
	print OUTFILE "$line\n";

}

#close files
close VIRAL;
close OUTFILE;
if ($bed) { close BEDFILE; }

#subroutines

sub printHelp {
	print "Pipeline for detection of viral integration sites within a genome\n\n";
	print "Usage:\n";
	print "\tperl softClip.pl --viral <sam> --thresh <f> --output <out>  --bed <bed> --help\n\n";
	print "Arguments:\n";
	print "\t--viral:   Alignment of reads to viral genomes (sam)\n";
	print "\t--thresh:  Threshold for fraction of read that must be aligned to not be vector rearrangement (0 < thresh < 1, default = 0.95)\n";
	print "\t--output:  Output file for results (default = rearrange.txt\n";
	print "\t--bed:		Output bed file (optional)\n";
	print "\t--verbose: Print progress messages\n";
	print "\t--help:    Print this help message and quit\n";

	exit;
}

