#!/usr/bin/perl
# It's not a bug, it's a feature

use strict;
use warnings;

use lib '.';

use ViralIntegration;
use Getopt::Long;

my $cutoff = 20; # default clipping cutoff
my $viral;
my $human;
my $output = "vecRearrange.txt";
my $verbose;

my $help;

GetOptions('cutoff=i' => \$cutoff,
		   'viral=s'  => \$viral,
		   'human=s'  => \$human,
		   'output=s' => \$output,
		   'verbose'  => \$verbose,
		   'help',    => \$help);

if ($help) { printHelp(); }

unless ($human and $viral) { printHelp(); }

my %viralReads;
my %humanReads;
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

	### Store informaton about the integration site:
	###	Integration Start Target = last clipped base position relative to target sequence
	###	Integration Stop Target  = first aligned base position relative to target sequence
	###	Integration Start Read 	 = last clipped base position relative to read length
	###	Integration Stop Read    = first aligned base position relative to read length
	###	Integration Orientation	 = orientation of integration site: + = after aligned sequence, - = before aligned sequence
	
	
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
	
	$viralReads{join("xxx", ($parts[0],$seq))} = join("xxx",($parts[2], $seq, $ori, $cig, $vSec, $vSup)); 

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

	my $cig = $parts[5]; # Process the CIGAR string to account for complex alignments
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

	if (exists $viralReads{join("xxx",($parts[0]), $seq)}) { # only consider reads that were tagged from the viral alignment, no need to consider excess reads
	
		my ($hSec, $hSup) = getSecSup($hl);

		$humanReads{join("xxx",($parts[0],$seq))} = join("\t",($parts[2], $seq, $ori, $cig, $hSec, $hSup))
	}
}
close HUMAN;

my ($rearrange, $line);
open (OUTFILE, ">$output") || die "Could not open output file: $output\n";
print OUTFILE "ReadID\tseq\tpossibleVecRearrange\n";


foreach my $key (keys %viralReads) {
	if (exists $humanReads{$key})	{
		#get info from alignments
		my $seq = (split("xxx",$key))[1];
		my ($otherpCig, $otherpDir) = (split("xxx",$humanReads{$key}))[3,2];
		my ($pCig, $pDir, $pSec, $pSup) = (split("xxx",$viralReads{$key}))[3,2,4,5];
		my $sup = join(";", $pSec, $pSup);
		my $readlen = length($seq);
		
		#get if possible vector rearrangement
		$line = join("\t", (split("xxx",$key))[0], (split("xxx",$key))[1], isRearrange($pCig, $pDir, $sup, $readlen, $otherpCig, $otherpDir, $seq));
		
		# print to file
		print OUTFILE "$line\n";
	}
}

close OUTFILE;

sub printHelp {
	print "Pipeline for detection of viral integration sites within a genome\n\n";
	print "Usage:\n";
	print "\tperl softClip.pl --viral <sam> --human <sam> --cutoff <n> --output <out> --bed <bed> --help\n\n";
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

