#!/usr/bin/perl
# It's not a bug, it's a feature

use strict;
use warnings;

use lib ".";

use ViralIntegration;
use Getopt::Long;

my $cutoff = 20; # default clipping cutoff
my $viral;
my $human;
my $output = "integrationSites.txt";
my $bed;
my $merged;
my $verbose;

my $help;

GetOptions('cutoff=i' => \$cutoff,
		   'viral=s'  => \$viral,
		   'human=s'  => \$human,
		   'output=s' => \$output,
		   'bed=s'    => \$bed,
		   'merged=s' => \$merged,
		   'verbose'  => \$verbose,
		   'help',    => \$help);

if ($help) { printHelp(); }

unless ($viral and $human) { printHelp(); }

### hash arrays will store putative discordant pairs from both alignments
### make one array for R1 and R2, and one for viral and human (four total)
my %viralR1;
my %viralR2;
my %humanR1;
my %humanR2;


### Collect junction reads from viral genome
### Junctions are defined by the CIGAR string: SM or MS
### Save the resulting information in a hash:
###	key: 	{ReadID}
###	value:	{sequence}xxx{TargetID}xxx{alignStart}xxx{SeqOrientation}xxx{CIGAR}
###		SeqOrientation: f = aligned in fwd orientation
###				r = aligned in rev orientation
###		Sequence is always given in fwd orientation if flag 0x10 is set, otherwise is reverse
### XA is the secondary alignments in the XA field (from BWA)
### SA is the supplementary alignments in the SA field (from BWA)

open (VIRAL, $viral) || die "Could not open viral alignment file: $viral\n";
if ($verbose) { print "Processing viral alignment...\n"; }
while (my $vl = <VIRAL>) {
	if ($vl =~ /^@/) { next; } # skip header lines

	my @parts = split(' ', $vl);
	
	# flags to consider
	#0x2 = read mapped in proper pair
	#0x4 = not mapped
	#0x8 = mate unmapped
	#0x10 = mapped in reverse orientation
	#0x40 = first in pair
	#0x80 = second in pair
	#0x800 = supplementary alignment

	if ($parts[1] & 0x800) { next(); } # skip supplementary alignments
	if ($parts[1] & 0x2) { next(); } # skip mapped in proper pair
	
	#want reads where if not mapped, mate is mapped, or vice versa
	
	unless ((($parts[1] & 0x8) and !($parts[1] & 0x4)) or (!($parts[1] & 0x8) and ($parts[1] & 0x4))) { next; }

	#get cigar
	my $cig = $parts[5]; # Don't process cigar for now; could be a cigar or * if unmapped
	
	#get sequence and orientation
	my ($seq, $seqOri);	
	if ($parts[1] & 0x10) { $seq = reverseComp($parts[9]); $seqOri = 'r'; }
	else { $seq = $parts[9]; $seqOri = 'f'; }
	
	#get read ID, aligned viral reference, alignment start
	my ($ID, $ref, $start) = ($parts[0], $parts[2], $parts[3]);
	
	#get supplementary (SA) and secondary (XA) alignments in order to check for possible vector rearrangements
	my ($vSec, $vSup) = getSecSup($vl);
	
	#get if first or second in pair
	my $inPair;
	if ($parts[1] & 0x40) 		{ $viralR1{$ID} = join("xxx", $seq, $seqOri, $ref, $start, $cig, $vSec, $vSup); }
	elsif ($parts[1] & 0x80) 	{ $viralR2{$ID} = join("xxx", $seq, $seqOri, $ref, $start, $cig, $vSec, $vSup); }

}
close VIRAL;

### Collect discordant reads from human genome
### This is the same process as for viral junctions
open (HUMAN, $human) || die "Could not open human alignment file: $human\n";
if ($verbose) { print "Processing human alignment...\n"; }
while (my $hl = <HUMAN>) {
	if ($hl =~ /^@/) { next; } # skip header lines

	my @parts = split(' ', $hl); #get fields from each alignment
	
	#want reads where if not mapped, mate is mapped, or vice versa
	unless ((($parts[1] & 0x8) and !($parts[1] & 0x4)) or (!($parts[1] & 0x8) and ($parts[1] & 0x4))) { next; }
	
	unless ((exists $viralR1{$parts[0]}) and (exists $viralR2{$parts[0]})) { next; } #only consider reads from viral alignment
	
	if ($parts[1] & 0x800) { next(); } # skip supplementary alignments
	if ($parts[1] & 0x2) { next(); } # skip mapped in proper pair
	
	#get cigar
	my $cig = $parts[5]; # Don't process cigar for now; could be a cigar or * if unmapped
	
	#get sequence and orientation
	my ($seq, $seqOri);	
	if ($parts[1] & 0x10) { $seq = reverseComp($parts[9]); $seqOri = 'r'; }
	else { $seq = $parts[9]; $seqOri = 'f'; }
	
	#get read ID, aligned viral reference, alignment start
	my ($ID, $ref, $start) = ($parts[0], $parts[2], $parts[3]);
	
	#get supplementary (SA) and secondary (XA) alignments in order to check for possible vector rearrangements
	my ($hSec, $hSup) = getSecSup($hl);
	
	#get if first or second in pair and add to appropriate array
	if ($parts[1] & 0x40) 		{ $humanR1{$ID} = join("xxx", $seq, $seqOri, $ref, $start, $cig, $hSec, $hSup);}
	elsif ($parts[1] & 0x80)  	{ $humanR2{$ID} = join("xxx", $seq, $seqOri, $ref, $start, $cig, $hSec, $hSup); }

}
close HUMAN;


### Look for evidence of integration sites by identifying discordant read-pairs
### Need to compare the junction sites in the viral and human genomes to see if there is any overlap (ambiguous bases)
### Integration coordinates are not exact, but assigned to the middle part of the human read
my @outLines;
if ($verbose) { print "Finding discordant read-pairs...\n"; }
foreach my $key (keys %viralR1) {

	#check that read ID is found for all hashes
	unless ( exists $viralR2{$key} ) { next; } 
	unless ( exists $humanR1{$key} ) { next; }
	unless ( exists $humanR2{$key} ) { next; }

	# Find discordant read pairs
	my @intData = findDiscordant($key, $viralR1{$key}, $viralR2{$key}, $humanR1{$key}, $humanR2{$key});
	
	# Once the positions are collected, put together the output	
	my $outLine;
	if ((scalar @intData) > 1) { $outLine = join("\t", @intData); }
	if (defined $outLine) { push(@outLines, join("\t", ($outLine, $key))) };
}

unless (@outLines) { die "No discordant read-pairs were detected\n"; } # if no integration events detected, finish

if ($verbose) { print "Writing output...\n"; }
printOutput($output, @outLines);
#if ($bed)    { printBed($bed, @outLines); }

#if ($merged) { printMerged($bed, $merged); }

exit;


#----------------------------------------------------------------------------------------------------------------------------------------
### Subroutines ###
#----------------------------------------------------------------------------------------------------------------------------------------

sub findDiscordant {
### Check if there is overlap between the integration sites
### in human and viral integration sites
### returns integration start/stop relative to read sequence

### BWA Alignments are 1-based
	my ($key, $vR1, $vR2, $hR1, $hR2) = @_;
	
	my ($seq1, $vR1ori, $vR1ref, $vR1start, $vR1cig, $vR1inPair, $vR1sec, $vR1sup) = (split('xxx', $vR1))[0, 1, 2, 3, 4, 5, 6, 7];
	my ($seq2, $vR2ori, $vR2ref, $vR2start, $vR2cig, $vR2inPair, $vR2sec, $vR2sup) = (split('xxx', $vR2))[0, 1, 2, 3, 4, 5, 6, 7];
	my ($hR1ori, $hR1ref, $hR1start, $hR1cig, $hR1inPair, $hR1sec, $hR1sup) = (split('xxx', $hR1))[1, 2, 3, 4, 5, 6, 7];
	my ($hR2ori, $hR2ref, $hR2start, $hR2cig, $hR2inPair, $hR2sec, $hR2sup) = (split('xxx', $hR2))[1, 2, 3, 4, 5, 6, 7];

	#check that we do have complementary discordant read-pairs in human and viral alignments
	unless ((( $vR1cig eq "*") and ($vR2cig =~ /^\d+M$/)) or (( $vR2cig eq "*") and ($vR1cig =~ /^\d+M$/))) { return; }
	unless ((( $hR1cig eq "*") and ($hR2cig =~ /^\d+M$/)) or (( $hR2cig eq "*") and ($hR1cig =~ /^\d+M$/))) { return; }
	if ((($vR1cig eq "*") and ($hR1cig eq "*")) or (($vR1cig =~ /^\d+M$/) and ($hR1cig =~ /^\d+M$/))) { return; }
	
	#get reference for virus and human based on which is mapped
	my ($hRef, $vRef, $hSeq, $vSeq);
	if ($hR1cig ne "*") { #if human matched R1
		$hRef = $hR1ref; 
		$hSeq = $seq1; 
		$vRef = $vR2ref;
		$vSeq = $seq2;
	}
	else { #if human matched R2
		$hRef = $hR2ref; 
		$hSeq = $seq2;
		$vRef = $vR1ref;
		$vSeq = $seq1; 
	}	
	
	my $readlen1 = length($seq1);
	my $readlen2 = length($seq2);
	
	#find if junction is human/virus or virus/human
	#find approximate location of junction on both human and virus side (assign to end of read closest to junction)
	#position indicated in SAM file is left-most mapping base: for forward reads this is 5' end, for reverse it's 3'
	my ($junct, $hIntStart, $hIntStop, $vIntStart, $vIntStop, $isVecRearrange, $isHumRearrange, $isHumAmbig, $isVirAmbig);
	if ($hR1cig ne "*") { #if R1 is matched in human
	
		#left-most mapped base is at position $hR1start, so need to add length of sequence to get position
		$hIntStart = $hR1start + $readlen1 - 1;
		$hIntStop = $hR1start + $readlen1;
		
		#opposite for viral read
		$vIntStart = $vR2start -1;
		$vIntStop = $vR2start;
		
		#get type of junction
		if ($hR1ori eq 'f') 	{ $junct = 'hv'; } #if human R1 is forward, orientation is hv
		elsif ($hR1ori eq 'r') 	{ $junct = 'vh'; } #otherwise orientation is vh
		
		#check for ambiguities
		#if ((join(";", $vR2sup, $vR2sec)) eq "NA;NA") { $isVecRearrange = "no"; }
		#else { $isVecRearrange = isRearrange($vR2cig, $vR2ori, (join(";", $vR2sup, $vR2sec)), $readlen2, $hR1cig, $hR1ori);}
	 
		#if ((join(";", $hR1sup, $hR1sec)) eq "NA;NA") { $isHumRearrange = "no"; }
		#else { $isHumRearrange = isRearrange($hR1cig, $hR1ori, (join(";", $hR1sup, $hR1sec)), $readlen1, $vR1cig, $vR1ori);}
	
		#check to see if location of human alignment is ambiguous: multiple equivalent alignments accounting for human part of read
		if ($hR1sec eq "NA") { $isHumAmbig = "no";}
		else { $isHumAmbig = isAmbigLoc($hR1ori, $hR1cig, $hR1sec);}
	
		#check to see if location of viral alignment is ambiguous: multiple equivalent alignments accounting for viral part of read
		if ($vR2sec eq "NA") { $isVirAmbig = "no";}
		else { $isVirAmbig = isAmbigLoc($vR2ori, $vR2cig, $vR2sec);}
	}
	elsif ($hR2cig ne '*') { #if R2 is matched in human
	
		#left-most mapped base is at postion $hR2start
		$hIntStart = $hR2start - 1;
		$hIntStop = $hR2start;
		
		#left-most mapped base is at position $vR1start, so need to add length of sequence to get position
		$vIntStart = $vR1start + $readlen1 -1;
		$vIntStop = $vR1start + $readlen1;
		
		#get type of junction
		if ($hR2ori eq 'r')		{ $junct = 'hv'; } #if human R2 is reverse, orientation is hv
		elsif ($hR2ori eq 'f') { $junct = 'vh'; } #if human R2 is forward, orientation is vh
	
		#check for ambiguities
		#if ((join(";", $vR1sup, $vR1sec)) eq "NA;NA") { $isVecRearrange = "no"; }
		#else { $isVecRearrange = isRearrange($vR1cig, $vR1ori, (join(";", $vR1sup, $vR1sec)), $readlen1, $hR2cig, $hR2ori);}
	 
		#if ((join(";", $hR2sup, $hR2sec)) eq "NA;NA") { $isHumRearrange = "no"; }
		#else { $isHumRearrange = isRearrange($hR2cig, $hR2ori, (join(";", $hR2sup, $hR2sec)), $readlen2, $vR2cig, $vR2ori);}
	
		#check to see if location of human alignment is ambiguous: multiple equivalent alignments accounting for human part of read
		if ($hR2sec eq "NA") { $isHumAmbig = "no";}
		else { $isHumAmbig = isAmbigLoc($hR2ori, $hR2cig, $hR2sec);}
	
		#check to see if location of viral alignment is ambiguous: multiple equivalent alignments accounting for viral part of read
		if ($vR1sec eq "NA") { $isVirAmbig = "no";}
		else { $isVirAmbig = isAmbigLoc($vR1ori, $vR1cig, $vR1sec);}
	}

	return($hRef, $hIntStart, $hIntStop, $vRef, $vIntStart, $vIntStop, '?', 'discordant', $junct, $hSeq, $vSeq, '-', join(';', $hR1sec, $hR2sec), join(';', $vR1sec, $vR2sec), 'NA', 'NA', $isVirAmbig, $isHumAmbig, 'NA');

}

sub printHelp {
	print "Pipeline for detection of viral integration sites within a genome\n\n";
	print "Usage:\n";
	print "\tperl discordant.pl --viral <sam> --human <sam> --cutoff <n> --output <out> --bed <bed> --help\n\n";
	print "Arguments:\n";
	print "\t--viral:   Alignment of reads to viral genomes (sam)\n";
	print "\t--human:   Alignment of reads to human genome (sam)\n";
	print "\t--cutoff:  Minimum number of clipped reads to be considered (default = 20)\n";
	print "\t--output:  Output file for results (default = integrationSite.txt\n";
	print "\t--bed:     Print integrations sites to indicated bed file (default = NA)\n";
	print "\t--merged:  Merge bedfile into overlapping integration sites (default = NA)\n";
	print "\t--verbose: Print progress messages\n";
	print "\t--help:    Print this help message and quit\n";

	exit;
}

