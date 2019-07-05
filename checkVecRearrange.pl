#!/usr/bin/perl
# It's not a bug, it's a feature

use strict;
use warnings;

use lib '.';

use ViralIntegration;
use Getopt::Long;

my $viral;
my $output = "rearrange.txt";
my $bed = "rearrange.bed";
my $verbose;
my $thresh = 0.95;
my $help;

GetOptions('thesh=f'  => \$thresh,
		   'viral=s'  => \$viral,
		   'output=s' => \$output,
		   'bed=s'	  => \$bed,
		   'verbose'  => \$verbose,
		   'help',    => \$help);

if ($help) { printHelp(); }

unless ($viral) { printHelp(); }

if (($thresh < 0) or ($thresh > 1)) { printHelp(); }

my %viralReads;
### These hashes contain the information about the integration sites within the resepctive genomes
### Data is saved as a tab seperated string
### Read name => Chr intStart intStop relStart relStop orientation readSeqeunce


### Collect junction reads from viral genome
### Junctions are defined by the CIGAR string: SM or MS
### Save the resulting information in a hash:
###	key: 	{ReadID}xxx{ReadSequence}
###	value:	{TargetID}xxx{Sequence}xxx{SeqOrientation}xxx{CIGAR}xxx{XA}xxx{SA}
###		SeqOrientation: f = aligned in fwd orientation
###				r = aligned in rev orientation
###		Sequence is always given in fwd orientation
### XA is the secondary alignments in the XA field (from BWA)
### SA is the supplementary alignments in the SA field (from BWA)

open (VIRAL, $viral) || die "Could not open viral alignment file: $viral\n";
if ($verbose) { print "Processing viral alignment...\n"; }
while (my $vl = <VIRAL>) {
	if ($vl =~ /^@/) { next; } # skip header lines

	my @parts = split("\t", $vl);
	

	if ($parts[1] & 0x800) { next(); } # skip supplementary alignments

	if ($parts[2] eq "*") { next; } # skip unaligned reads

	my $cig = $parts[5]; 
	unless ($cig) { next; } # keep checking to make sure double clipped reads don't sneak through

	#get information about location of alignment
	my $ref = $parts[2];
	my $pos = $parts[3];
	
	
	#get supplementary (SA) and secondary (XA) alignments in order to check for possible vector rearrangements
	my ($vSec, $vSup) = getSecSup($vl);
	
	my $seq;
	my $ori;
	if ($parts[1] & 0x10) { 
		($seq) = reverseComp($parts[9]); 
		$ori   = 'r';
	}
	else				  { 
		$seq = $parts[9]; 
		$ori = 'f';
	}
	
	$viralReads{join("xxx", ($parts[0],$seq))} = join("xxx",($parts[2], $seq, $ori, $cig, $vSec, $vSup, $pos)); 

}
close VIRAL;

my ($rearrange, $line, $bedLine);
open (OUTFILE, ">$output") || die "Could not open output file: $output\n";
if ($bed) { open (BEDFILE, ">$bed") || die "Could not open output file: $output\n"; }
print OUTFILE "ReadID\tpossibleVecRearrange\tTotalGapBP\tGaps\tsups+secs\tseq\n";


foreach my $key (keys %viralReads) {
		#get info from alignments
		my $seq = (split("xxx",$key))[1];
		my ($pCig, $pDir, $pSec, $pSup, $pos, $ref) = (split("xxx",$viralReads{$key}))[3,2,4,5,6,0];
		my $sup = join(";", $pSec, $pSup);
		
		#get if possible vector rearrangement
		my ($isRearrange, $gapBP, $gaps) = isRearrange($pCig, $pDir, $ref, $pos, $sup, $seq, $thresh);
		$line = join("\t", (split("xxx",$key))[0], $isRearrange, $gapBP, $gaps, $sup, $seq);
		
		print OUTFILE "$line\n";
		
		# print to file
		if ($bed) {
			#first get info about primary alignment
			my @aligns = getMatchedRegions($pCig, $pDir);
			
			foreach my $align (@aligns) {
				
			
			}
		
			#print BEDFILE "$bedLine\n";
		}
		
		
		
		#print to bed file
		
}

close OUTFILE;

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

