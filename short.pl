#!/usr/bin/perl
# It's not a bug, it's a feature

#this script identifies short insertions of viral DNA into a host genome
#making use of an alignment (SAM file) of NGS data against both the host and viral genomes

#identity insertions as reads which contain soft-clips on each side in the viral alignment
#and contain matched regions on each side and insertion in the middle in host alignment

use strict;
use warnings;

use lib '.';

use ViralIntegration;
use Getopt::Long;

my $cutoff = 20; # default clipping cutoff
my $thresh = 0.95; #default amount of read that must be covered by alignments for rearrangement
my $tol = 5; #when processing CIGARS, combine any IDPN elements between matched regions with this number of bases or less
my $viral;
my $human;
my $output = "integrationSites.txt";
my $bed;
my $merged;
my $verbose;

my $help;

GetOptions('cutoff=i' => \$cutoff,
		   'thresh=f' => \$thresh,
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
	if ($parts[1] & 0x4) {next; } #skip unmapped alignments

	if ($parts[2] eq "*") { next; } # skip unaligned reads

	my $cig = processCIGAR2($parts[5], $tol); # Process the CIGAR string to account for small insertions/deletions
	unless ($cig) { next; } # keep checking to make sure double clipped reads don't sneak through
	
	unless ($cig =~ /(^\d+[SH].*\d+[SH]$)/) { next; } # want alignments where both ends are clipped
	
	#number of mapped bases must be more than $cutoff
	my (@insert) = ($cig =~ /(\d+)M/g);
	my $mapBP = eval join("+", @insert);
	if ($mapBP < $cutoff) { next; }
	
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

	#get read ID, aligned viral reference, alignment start
	my ($ID, $ref, $start) = ($parts[0], $parts[2], $parts[3]);
	
	#get supplementary (SA) and secondary (XA) alignments in order to check for possible vector rearrangements
	my ($vSec, $vSup) = getSecSup($vl);
	
	$viralIntegrations{join("xxx",($parts[0],$seq))} = join("xxx", $seq, $ori, $ref, $start, $cig, $vSec, $vSup); 
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
	if ($parts[1] & 0x4) {next; } #skip unmapped alignments
	
	unless ($parts[5]) { next; }

	my $cig = processCIGAR2($parts[5], $tol); # Process the CIGAR string to account for complex alignments
	unless ($cig) { next; }
	
	unless ($cig =~ /(^\d+M.+\d+M$)/) { next; } #alignments must be matched at both ends
	unless ($cig =~ /I/) { next; } #must have inserted region
	
	#only want alignments where matched region is more than cutoff
	my ($match1, $match2) = ($cig =~ /^(\d+)M.+[ISDNPH](\d+)M$/);
	if (($match1 < $cutoff) or ($match2 < $cutoff)) { next; }
	
	#also want alignments where inserted region is more than cutoff
	my (@insert) = ($cig =~ /(\d+)I/g);

	
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
	
	#get read ID, aligned viral reference, alignment pos
	my ($ID, $ref, $start) = ($parts[0], $parts[2], $parts[3]);

	if (exists $viralIntegrations{join("xxx",($parts[0],$seq))}) { # only consider reads that were tagged from the viral alignment, no need to consider excess reads
	
		my ($hSec, $hSup) = getSecSup($hl);
		
		$humanIntegrations{join("xxx",($parts[0],$seq))} = join("xxx", $seq, $ori, $ref, $start, $cig, $hSec, $hSup);
	}
}
close HUMAN;


### Look for evidence of integration sites by identifying chimeric reads
### Need to compare the junction sites in the viral and human genomes to see if there is any overlap (ambiguous bases)
### Adjust integration coordinates appropriately
my $count = 0;
my @outLines;
if ($verbose) { print "Detecting chimeric reads...\n"; }
foreach my $key (keys %viralIntegrations) {
	if (exists $humanIntegrations{$key}) { # only consider reads that are flagged in both human and viral alignments
	
		$count += 1;
		
		#as a first step, just save any potential short integrations in a txt file
		my @viralData = split('xxx', $viralIntegrations{$key});
		my @hostData = split('xxx', $humanIntegrations{$key});
		my ($ID, $seq) = split('xxx', $key);
		
		push(@outLines, join("\t", $hostData[2], $hostData[3], $hostData[4], $viralData[2], $viralData[3], $viralData[4], $ID, $seq));
		
	}
}


#print file
my $header = join("\t", "HostRef", "HostPos", "HostCig", "ViralRef", "ViralPos", "ViralCig", "ID", "seq\n");

if ($verbose) { print "Writing output...\n"; }

printOutput($output, $header, @outLines); #write to outfile: if no sites detected will be header only


exit;

#----------------------------------------------------------------------------------------------------------------------------------------
### Subroutines ###
#----------------------------------------------------------------------------------------------------------------------------------------

sub printHelp {
	print "Pipeline for detection of viral integration sites within a genome\n\n";
	print "Usage:\n";
	print "\tperl short.pl --viral <sam> --human <sam> --cutoff <n> --thresh <n> --output <out> --bed <bed> --help\n\n";
	print "Arguments:\n";
	print "\t--viral:   Alignment of reads to viral genomes (sam)\n";
	print "\t--human:   Alignment of reads to human genome (sam)\n";
	print "\t--cutoff:  Minimum number of clipped reads to be considered (default = 20)\n";
	print "\t--thresh:	Amount of read that must be covered by one alignment to be considered rearrangement";
	print "\t--output:  Output file for results (default = integrationSite.txt\n";
	print "\t--bed:     Print integrations sites to indicated bed file (default = NA)\n";
	print "\t--merged:  Merge bedfile into overlapping integration sites (default = NA)\n";
	print "\t--verbose: Print progress messages\n";
	print "\t--help:    Print this help message and quit\n";

	exit;
}

