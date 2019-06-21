#!/usr/bin/perl
# It's not a bug, it's a feature

use strict;
use warnings;

use File::Basename qw(dirname);
use Cwd  qw(abs_path);
use lib dirname(dirname abs_path $0) . '/lib';

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

my %viralIntegrations;
my %humanIntegrations;
### These hashes contain the information about the integration sites within the resepctive genomes
### Data is saved as a tab seperated string
### Read name => Chr intStart intStop relStart relStop orientation readSeqeunce


### Collect junction reads from viral genome
### Junctions are defined by the CIGAR string: SM or MS
### Save the resulting information in a hash:
###	key: 	{ReadID}xxx{ReadSequence}
###	value:	{TargetID}xxx{ViralIntegrationData}xxx{Sequence}xxx{SeqOrientation}xxx{CIGAR}xxx{XA}xxx{SA}
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
	if ($parts[5] =~ /(^\d+[SH].*\d+[SH]$)/) { next; } # skip alignments where both ends are clipped

	unless ($parts[5] =~ /^\d+[SH]|\d+[SH]$/) { next; } # at least one end of the read must be clipped in order to be considered

	my $cig = processCIGAR($parts[5], $parts[9]); # Process the CIGAR string to account for complex alignments
	unless ($cig) { next; } # keep checking to make sure double clipped reads don't sneak through

	### Store informaton about the integration site:
	###	Integration Start Target = last clipped base position relative to target sequence
	###	Integration Stop Target  = first aligned base position relative to target sequence
	###	Integration Start Read 	 = last clipped base position relative to read length
	###	Integration Stop Read    = first aligned base position relative to read length
	###	Integration Orientation	 = orientation of integration site: + = after aligned sequence, - = before aligned sequence
	my @viralInt;
	if    ($cig =~ /^\d+[SH].+\d+[SH]$/) { @viralInt = undef; } # double check, sometimes double clipped sneak through
	elsif ($cig =~ /(^\d+)[SH]/)   	  { if ($1 > $cutoff) { @viralInt = analyzeRead($parts[3], $cig, "-"); } } # integration site is before the viral sequence
	elsif ($cig =~ /(\d+)[SH]$/) 	  { if ($1 > $cutoff) { @viralInt = analyzeRead($parts[3], $cig, "+"); } } # integration site is after the viral sequence

	# ID = readName_seq (in forward direction)

	#get supplementary (SA) and secondary (XA) alignments in order to check for possible vector rearrangements
	
	my ($vSec, $vSup) = getSecSup($vl);
	
	
	if   (@viralInt and ($parts[1] & 0x10)) { $viralIntegrations{join("xxx", ($parts[0],(reverseComp($parts[9]))[0]))} = join("\t",($parts[2], @viralInt, (reverseComp($parts[9]))[0], 'r', $cig, $vSec, $vSup)); }
	elsif (@viralInt) 			{ $viralIntegrations{join("xxx", ($parts[0],$parts[9]))}	           = join("\t",($parts[2], @viralInt, $parts[9],                   'f', $cig, $vSec, $vSup)); } 
	#$viralIntegrations{join("xxx", ($parts[0],$parts[9]))} = join("\t",($parts[2], @viralInt, $parts[9]));
}
close VIRAL;

### Collect junction reads from human genome
### This is the same process as for viral junctions
open (HUMAN, $human) || die "Could not open human alignment file: $human\n";
if ($verbose) { print "Processing human alignment...\n"; }
while (my $hl = <HUMAN>) {
	if ($hl =~ /^@/) { next; } # skip header lines

	my @parts = split("\t", $hl); #get fields from each alignment
	
	if ($parts[1] & 0x800) { next(); } # skip supplementary alignments

	if ($parts[2] eq "*") { next; } # skip unaligned reads
	unless ($parts[5]) { 
		next; 
	}
	if ($parts[5] =~ /(^\d+[SH].*\d+[SH]$)/) { next; } #skip alignments where both ends are clipped
	
	unless ($parts[5] =~ /^\d+[SH]|\d+[SH]$/) { next; }
	my $cig = processCIGAR($parts[5], $parts[9]); # Process the CIGAR string to account for complex alignments
	unless ($cig) { next; }

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

	#$seq = $parts[9];

	if (exists $viralIntegrations{join("xxx",($parts[0],$seq))}) { # only consider reads that were tagged from the viral alignment, no need to consider excess reads
	
		my ($hSec, $hSup) = getSecSup($hl);
		
		my @humanInt;
		if    ($cig =~ /^\d+[SH].+\d+[SH]$/) { @humanInt = undef; }
		elsif ($cig =~ /(^\d+)[SH]/)         { if ($1 > $cutoff) { @humanInt = analyzeRead($parts[3], $cig, "-"); } }
		elsif ($cig =~ /(\d+)[SH]$/)         { if ($1 > $cutoff) { @humanInt = analyzeRead($parts[3], $cig, "+"); } }

		if (@humanInt) { $humanIntegrations{join("xxx",($parts[0],$seq))} = join("\t",($parts[2], @humanInt, $seq, $ori, $cig, $hSec, $hSup)); }
	}
}
close HUMAN;


### Look for evidence of integration sites by identifying discordant read-pairs
### Need to compare the junction sites in the viral and human genomes to see if there is any overlap (ambiguous bases)
### Integration coordinates are not exact, but assigned to the middle part of the human read
my @outLines;
if ($verbose) { print "Detecting integrations...\n"; }
foreach my $key (keys %viralIntegrations) {
	if (exists $humanIntegrations{$key}) { # only consider reads that are flagged in both human and viral alignments
		# Collect positions relative to read
		# Returns an array with the following:
		#	Integration Start
		#	Integration Stop
		#	Human Start
		#	Human Stop
		#	Viral Start
		#	Viral Stop	
		#	Overlap
		my @intData = findDiscordant($viralIntegrations{$key}, $humanIntegrations{$key});
		# Once the positions are collected, put together the output	
		my $outLine;
		if (@intData) { $outLine = extractOutput($viralIntegrations{$key}, $humanIntegrations{$key}, @intData); }
		if ($outLine) { push(@outLines, join("\t", ($outLine, (split("xxx",$key))[0],(split("xxx",$key))[1]))) };
	}
}

unless (@outLines) { die "No integration sites were detected\n"; } # if no integration events detected, finish

if ($verbose) { print "Writing output...\n"; }
printOutput($output, @outLines);
if ($bed)    { printBed($bed, @outLines); }

if ($merged) { printMerged($bed, $merged); }

exit;


#----------------------------------------------------------------------------------------------------------------------------------------
### Subroutines ###
#----------------------------------------------------------------------------------------------------------------------------------------


sub findDiscordant {
### Check if there is overlap between the integration sites
### in human and viral integration sites
### returns integration start/stop relative to read sequence

### BWA Alignments are 1-based
	my ($viralData, $humanData) = @_;

	my ($vStart, $vStop, $vOri, $seq, $vDir, $vCig, $vSec, $vSup) = (split("\t",$viralData))[3,4,5,6,7,8,-2,-1];
	my ($hStart, $hStop, $hOri, $hDir, $hCig, $hSec, $hSup)       = (split("\t",$humanData))[3,4,5,7,8,-2,-1];

	if	((($vOri eq $hOri) and ($vDir eq $hDir)) or( ($vOri ne $hOri) and ($vDir ne $hDir))) { return; } # in some cases the same part of the read may be clippeed 
	### CIGAR strings are always reported relative to the strand
	### 100M50S on a fwd read = 50S100M on a rev read
	### Therefore if integration position and read orientation match (e.g. ++ and ff) same part of the read is clipped
	### If they all don't match (e.g. +- and fr which is equivalent to ++ and ff) then the same part of the read is clipped
	### ++ and fr is ok because this is the same as +- and ff

	unless 	($vOri and $hOri) { return; } # Catch a weird case where an orientation isn't found. Appears to happen when both ends are clipped

	### First find if there is any overlap between the human and viral junctions
	### Do so by comparing the CIGAR strings
	### Can calculate by subtracting aligned of one from clipped of other	
	my ($hClip) = ($hCig =~ /(\d+)S/);
	my ($hAlig) = ($hCig =~ /(\d+)M/);
	my ($vClip) = ($vCig =~ /(\d+)S/);
	my ($vAlig) = ($vCig =~ /(\d+)M/);
	
	### Overlap should be the same regardless of how it's calculated so double check
	my $overlap1 = abs($hAlig - $vClip);
	my $overlap2 = abs($vAlig - $hClip);

	#here check that absolute values of overlap are the same
	unless (abs($overlap1) == abs($overlap2)) { 
		print "Impossible overlap found\n";
		return;
	}
	
	#ambigous bases may result from either an overlap of aligned regions, or a gap
	#if the total number of aligned bases (from human and viral alignments) are greater than the read length, then it's an overlap
	#otherwise it's a gap
	my ($overlaptype, $readlen);
	$readlen = length($seq);
	if 		(($hAlig + $vAlig) > ($readlen))  {	$overlaptype = "overlap";	} #overlap
	elsif 	(($hAlig + $vAlig) == ($readlen)) {	$overlaptype = "none";		} #no gap or overlap
	else 									  {	$overlaptype = "gap";		} #gap
	
	
	#check to see is whole read can be accounted for by vector rearrangements
	
	my $isVecRearrange;
	if ((join(";", $vSup, $vSec)) eq "NA;NA") { $isVecRearrange = "no"; }
	else { $isVecRearrange = isRearrange($vCig, $vDir, (join(";", $vSup, $vSec)), $readlen);}
	 
	my $isHumRearrange;
	if ((join(";", $hSup, $hSec)) eq "NA;NA") { $isHumRearrange = "no"; }
	else { $isHumRearrange = isRearrange($hCig, $hDir, (join(";", $hSup, $hSec)), $readlen);}
	
	
	#check to see if location of human alignment is ambiguous: multiple equivalent alignments accounting for human part of read
	my $isHumAmbig;
	if ($hSec eq "NA") { $isHumAmbig = "no";}
	else { $isHumAmbig = isAmbigLoc($hDir, $hCig, $hSec, $hAlig);}
	
	#check to see if location of viral alignment is ambiguous: multiple equivalent alignments accounting for viral part of read
	my $isVirAmbig;
	if ($vSec eq "NA") { $isVirAmbig = "no";}
	else { $isVirAmbig = isAmbigLoc($vDir, $vCig, $vSec, $vAlig);}

	### Calculate the start and stop positions of the viral and human sequences relative to the read
	my ($hRStart,$hRStop) = extractSeqCoords($hOri, $hDir, $hAlig, abs($overlap1), $readlen);
	my ($vRStart,$vRStop) = extractSeqCoords($vOri, $vDir, $vAlig, abs($overlap1), $readlen);

	### Collect integration start/stop positions relative to read
	### These are the bases that flank the bond that is broken by insertion of the virus
	### Takes any overlap into account	
	my ($intRStart, $intRStop, $order);

	if    ($hRStop < $vRStart) { ($intRStart,$intRStop,$order) = ($hRStop,$vRStart,"hv"); } # Orientation is Human -> Virus
	elsif ($vRStop < $hRStart) { ($intRStart,$intRStop,$order) = ($vRStop,$hRStart,"vh"); } # Orientation is Virus -> Human
	else 			   { 
		print "Something weird has happened";
		return;
	}
	

	return($intRStart, $intRStop, $hRStart, $hRStop, $vRStart, $vRStop, $overlap1, $order, $overlaptype, $isVecRearrange, $isHumRearrange, $isHumAmbig, $isVirAmbig);

	#if 	  ($vOri eq "+" and $hOri eq "-" and $vStop - 1 > $hStart) { return($hStart, $vStop); } # overalp with order virus -> human
	#elsif ($vOri eq "-" and $hOri eq "+" and $hStop - 1 > $vStart)    { return($vStart, $hStop); } # overlap with order human -> virus
	#else 							   	   { return($hStart, $hStop); } # if no overlap, use human positions
}

sub printHelp {
	print "Pipeline for detection of viral integration sites within a genome\n\n";
	print "Usage:\n";
	print "\tperl detectionPipeline.pl --viral <sam> --human <sam> --cutoff <n> --output <out> --bed <bed> --help\n\n";
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