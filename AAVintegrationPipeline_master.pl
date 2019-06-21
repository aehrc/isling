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

### Look for valid integration sites
### Need to compare the junction sites in the viral and human genomes to see if there is any overlap (ambiguous bases)
### Adjust integration coordinates appropriately
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
		my @intData = collectIntersect($viralIntegrations{$key}, $humanIntegrations{$key});
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

