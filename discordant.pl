#!/usr/bin/perl
# It's not a bug, it's a feature

use strict;
use warnings;

use lib ".";

use ViralIntegration;
use Getopt::Long;

my $cutoff = 20; # default clipping cutoff
my $thresh = 0.95; #default amount of read that must be covered by alignments for rearrangement
my $tol = 3; #when processing CIGARS, combine any IDPN elements between matched regions with this number of bases or less
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
	
	#for read pairs where one read is mapped and the other is not:
	#want reads where if not mapped, mate is mapped, or vice versa
	
	#but if allowing some clipping, might have case were one read has small portion mapped facing into the junction
	#so in this case also need to consider reads where both mates are mapped
	
	#so essentially only exclude reads where both reads are unmapped
	
	if (($parts[1] & 0x8) and ($parts[1] & 0x4)) { next; }

	#get cigar
	my ($cig, $editDist2) = processCIGAR2($parts[5], $tol); # Note that could be a cigar or * if unmapped
	
	
	#get sequence and orientation
	my ($seq, $seqOri);	
	if ($parts[1] & 0x10) { $seq = reverseComp($parts[9]); $seqOri = 'r'; }
	else { $seq = $parts[9]; $seqOri = 'f'; }
	
	#get read ID, aligned viral reference, alignment start
	my ($ID, $ref, $start) = ($parts[0], $parts[2], $parts[3]);
	
	#get supplementary (SA) and secondary (XA) alignments in order to check for possible vector rearrangements
	my ($vSec, $vSup) = getSecSup($vl);
	my $editDist = getEditDist($vl);

	#if the read is not matched, cant get edit distance  - just set it to 0
	unless ($editDist) { $editDist = 0; $editDist2 = 0; } 
	
	#get if first or second in pair
	if ($parts[1] & 0x40) 		{ $viralR1{$ID} = join("xxx", $seq, $seqOri, $ref, $start, $cig, $vSec, $vSup, ($editDist+$editDist2)); }
	elsif ($parts[1] & 0x80) 	{ $viralR2{$ID} = join("xxx", $seq, $seqOri, $ref, $start, $cig, $vSec, $vSup, ($editDist+$editDist2)); }

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
	if (!($parts[1] & 0x8) and !($parts[1] & 0x4)) { next; }
	
	unless ((exists $viralR1{$parts[0]}) and (exists $viralR2{$parts[0]})) { next; } #only consider reads from viral alignment
	
	if ($parts[1] & 0x800) { next(); } # skip supplementary alignments
	if ($parts[1] & 0x2) { next(); } # skip mapped in proper pair
	
	#get cigar
	my ($cig, $editDist2) = processCIGAR2($parts[5], $tol); # Note that could be a cigar or * if unmapped
	
	#get sequence and orientation
	my ($seq, $seqOri);	
	if ($parts[1] & 0x10) { $seq = reverseComp($parts[9]); $seqOri = 'r'; }
	else { $seq = $parts[9]; $seqOri = 'f'; }
	
	#get read ID, aligned viral reference, alignment start
	my ($ID, $ref, $start) = ($parts[0], $parts[2], $parts[3]);
	
	#get supplementary (SA) and secondary (XA) alignments in order to check for possible vector rearrangements
	my ($hSec, $hSup) = getSecSup($hl);
	my $editDist = getEditDist($hl);
	
	#if the read is not matched, cant get edit distance  - just set it to 0
	unless ($editDist) { $editDist = 0; $editDist2 = 0; } 
	
	#get if first or second in pair and add to appropriate array
	if ($parts[1] & 0x40) 		{ $humanR1{$ID} = join("xxx", $seq, $seqOri, $ref, $start, $cig, $hSec, $hSup, ($editDist+$editDist2));}
	elsif ($parts[1] & 0x80)  	{ $humanR2{$ID} = join("xxx", $seq, $seqOri, $ref, $start, $cig, $hSec, $hSup, ($editDist+$editDist2)); }

}
close HUMAN;


### Look for evidence of integration sites by identifying discordant read-pairs
### Need to compare the junction sites in the viral and human genomes to see if there is any overlap (ambiguous bases)
### Integration coordinates are not exact, but assigned to the middle part of the human read
my @outLines;
if ($verbose) { print "Finding discordant read-pairs...\n"; }
foreach my $key (keys %humanR1) {

	#check that read ID is found for all hashes
	unless ( exists $viralR1{$key} ) { next; } 
	unless ( exists $viralR2{$key} ) { next; }
	unless ( exists $humanR2{$key} ) { next; }

	# Find discordant read pairs
	my @intData = findDiscordant($key, $viralR1{$key}, $viralR2{$key}, $humanR1{$key}, $humanR2{$key});
	unless (@intData) { next; }	

	# Once the positions are collected, put together the output	
	my $outLine = join("\t", @intData);
	my $combReads = join("xxx", (split("xxx", $viralR1{$key}))[0], (split("xxx", $viralR2{$key}))[0]);
	if (defined $outLine) { push(@outLines, join("\t", ($outLine, $key, $combReads))); }
}



if ($verbose) { print "Writing output...\n"; }

my $header = "Chr\tIntStart\tIntStop\tVirusRef\tVirusStart\tVirusStop\tNoAmbiguousBases\tOverlapType\tOrientation\tHostSeq\tViralSeq\tAmbiguousSeq\tHostEditDist\tViralEditDist\tTotalEditDist\tPossibleHostTranslocation\tPossibleVectorRearrangement\tHostPossibleAmbiguous\tViralPossibleAmbiguous\tReadID\tmerged\n";		
		
printOutput($output, $header, @outLines); #write to outfile: if no sites detected will be header only

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
	my ($key, $vR1, $vR2, $hR1, $hR2) = @_;
	
	my ($seq1, $vR1ori, $vR1ref, $vR1start, $vR1cig, $vR1sec, $vR1sup, $vNM1) = (split('xxx', $vR1))[0, 1, 2, 3, 4, 5, 6, -1];
	my ($seq2, $vR2ori, $vR2ref, $vR2start, $vR2cig, $vR2sec, $vR2sup, $vNM2) = (split('xxx', $vR2))[0, 1, 2, 3, 4, 5, 6, -1];
	my ($hR1ori, $hR1ref, $hR1start, $hR1cig, $hR1sec, $hR1sup, $hNM1) = (split('xxx', $hR1))[1, 2, 3, 4, 5, 6, -1];
	my ($hR2ori, $hR2ref, $hR2start, $hR2cig, $hR2sec, $hR2sup, $hNM2) = (split('xxx', $hR2))[1, 2, 3, 4, 5, 6, -1];

	#discordant read-pair can be either one mate mapped and one unmapped
	#or one mate mostly mapped and the other with only a small mapped region
	
	#find if R1 and R2 are mostly mapped or unmapped
	my ($vR1map, $vR1mapBP) = isMapped($vR1cig, $cutoff);
	my ($vR2map, $vR2mapBP) = isMapped($vR2cig, $cutoff);
	my ($hR1map, $hR1mapBP) = isMapped($hR1cig, $cutoff);
	my ($hR2map, $hR2mapBP) = isMapped($hR2cig, $cutoff);
	
	#check that we do have complementary discordant read-pairs in human and viral alignments
	#check that both reads are not either both mapped or unmapped in each alignment
	unless (($vR1map ne $vR2map) and ($hR1map ne $hR2map)) { return; }
	unless (($vR1map eq $hR2map) and ($hR2map eq $vR1map)) { return; }
	
	#if a read-pair has gotten to this point, check that number of unmapped bases in 'map'
	#is less than $cutoff
	#do this to make sure we're reasonably confident of mapping and avoid reads with relatively
	#small mapped regions (eg 91S25M34S) getting through
	
	my $readlen1 = length($seq1);
	my $readlen2 = length($seq2);
	
	#get reference for virus and human based on which is mapped
	my ($hRef, $vRef, $hSeq, $vSeq, $hNM, $vNM);
	if ($hR1map eq "map") { #if human matched R1
		unless (($readlen1 - $hR1mapBP) < $cutoff) { return; } #check unmapped bases is less than cutoff
		unless (($readlen2 - $vR2mapBP) < $cutoff) { return; } #check unmapped bases is less than cutoff
		$hRef = $hR1ref; 
		$hSeq = $seq1; 
		$hNM = $hNM1;
		$vRef = $vR2ref;
		$vSeq = $seq2;
		$vNM = $vNM2;
	}
	else { #if human matched R2
		unless (($readlen2 - $hR2mapBP) < $cutoff) { return; } #check unmapped bases is less than cutoff
		unless (($readlen1 - $vR1mapBP) < $cutoff) { return; } #check unmapped bases is less than cutoff
		$hRef = $hR2ref; 
		$hSeq = $seq2;
		$hNM = $hNM2;
		$vRef = $vR1ref;
		$vSeq = $seq1; 
		$vNM = $vNM1;
	}	
	
	#find if junction is human/virus or virus/human
	#find approximate location of junction on both human and virus side (assign to end of read closest to junction)
	#position indicated in SAM file is left-most mapping base: for forward reads this is 5' end, for reverse it's 3'
	my ($junct, $hIntStart, $hIntStop, $vIntStart, $vIntStop, $isVecRearrange, $isHumRearrange, $isHumAmbig, $isVirAmbig);
	if ($hR1map eq "map") { #if R1 is matched in human
		
		#get integration start and stop positions
		($hIntStart, $hIntStop, $vIntStart, $vIntStop) = getIntPos($hR1cig, $hR1ori, $hR1start, $readlen1, $vR2cig, $vR2ori, $vR2start, $readlen2);

		#get type of junction
		if ($hR1ori eq 'f') 	{ $junct = 'hv'; } #if human R1 is forward, orientation is hv
		elsif ($hR1ori eq 'r') 	{ $junct = 'vh'; } #otherwise orientation is vh
		
		#check for ambiguities
		if ((join(";", $vR2sup, $vR2sec)) eq "NA;NA") { $isVecRearrange = "no"; }
		else { $isVecRearrange = (isRearrange($vR2cig, $vR2ori, $vR2ref, $vR2start, (join(";", $vR2sup, $vR2sec)), $seq2, $thresh))[0];}
		
		if ((join(";", $hR1sec, $hR1sup)) eq "NA;NA") { $isHumRearrange = "no"; }
		else { $isHumRearrange = (isRearrange($hR1cig, $hR1ori, $hR1ref, $hR1start, (join(";", $hR1sec, $hR1sup)), $seq1, $thresh))[0];}
	 
		
	
		#check to see if location of human alignment is ambiguous: multiple equivalent alignments accounting for human part of read
		if ($hR1sec eq "NA") { $isHumAmbig = "no";}
		else { $isHumAmbig = isAmbigLoc($hR1ori, $hR1cig, $hR1sec, 'discordant', $seq1, $tol);}
	
		#check to see if location of viral alignment is ambiguous: multiple equivalent alignments accounting for viral part of read
		if ($vR2sec eq "NA") { $isVirAmbig = "no";}
		else { $isVirAmbig = isAmbigLoc($vR2ori, $vR2cig, $vR2sec, 'discordant', $seq2, $tol);}
	}
	elsif ($hR2map eq 'map') { #if R2 is matched in human
	
		#get integration start and stop positions
		($vIntStart, $vIntStop, $hIntStart, $hIntStop) = getIntPos($vR1cig, $vR1ori, $vR1start, $readlen1, $hR2cig, $hR2ori, $hR2start, $readlen2);
		
		#get type of junction
		if ($hR2ori eq 'r')		{ $junct = 'hv'; } #if human R2 is reverse, orientation is hv
		elsif ($hR2ori eq 'f') { $junct = 'vh'; } #if human R2 is forward, orientation is vh
	
		#check for ambiguities
		if ((join(";", $vR1sup, $vR1sec)) eq "NA;NA") { $isVecRearrange = "no"; }
		else { $isVecRearrange = (isRearrange($vR1cig, $vR1ori, $vR1ref, $vR1start, (join(";", $vR1sup, $vR1sec)), $seq1, $thresh))[0];}
		
		if ((join(";", $hR2sup, $hR2sec)) eq "NA;NA") { $isHumRearrange = "no"; }
		else { $isHumRearrange = (isRearrange($hR2cig, $hR2ori, $hR2ref, $hR2start, (join(";", $hR2sup, $hR2sec)), $seq2, $thresh))[0];}
	 	
		#check to see if location of human alignment is ambiguous: multiple equivalent alignments accounting for human part of read
		if ($hR2sec eq "NA") { $isHumAmbig = "no";}
		else { $isHumAmbig = isAmbigLoc($hR2ori, $hR2cig, $hR2sec, 'discordant', $seq2, $tol);}
	
		#check to see if location of viral alignment is ambiguous: multiple equivalent alignments accounting for viral part of read
		if ($vR1sec eq "NA") { $isVirAmbig = "no";}
		else { $isVirAmbig = isAmbigLoc($vR1ori, $vR1cig, $vR1sec, 'discordant', $seq1, $tol);}
	}
		return($hRef, $hIntStart, $hIntStop, $vRef, $vIntStart, $vIntStop, '?', 'discordant', $junct, $hSeq, $vSeq, '-', $hNM, $vNM, $isHumRearrange, $isVecRearrange, $isVirAmbig, $isHumAmbig);

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

sub isMapped{

	my ($cig, $cutoff) = @_;
	
	#check if cigar is 'mapped' or not
	
	#read is considered to be unmapped if number of mapped bases is less than $cutoff
	#otherwise it's mapped
	
	#return 'map' or 'unmap' for mapped or unmapped read respectively
	#also return number of mapped bases
	
	if ($cig =~ /\*/) { return ('unmap', 0); }
	
	my (@letters, @numbers);
	getCigarParts($cig, \@letters, \@numbers);
	
	my $mappedBP =0 ;
	for my $i (0..$#letters) {
		if ($letters[$i] =~ /M/) { $mappedBP += $numbers[$i]; }
	}

	if ($mappedBP > $cutoff) { return ('map', $mappedBP); }
	return ('unmap', $mappedBP);
	
}

sub getIntPos {

	my ($mapR1cig, $mapR1ori, $mapR1start, $readlen1, $mapR2cig, $mapR2ori, $mapR2start, $readlen2) = @_;
	#get positions of last mapped base facing in towards junction
	#
	
	#regardless of whether R1 is mapped in reverse or forward orientation,
	#the estimate of int site position will always be pos (1-based left-most mapping postion) 
	#plus the read length
	my $R1IntStart = $mapR1start + getLastMapped($mapR1cig, $mapR1ori) - 1;
	my $R1IntStop = $R1IntStart;

	#similarly, estimate for R2 is always at pos (1-based mapping position)
	my $R2IntStart = $mapR2start;
	my $R2IntStop = $R2IntStart;

	return ($R1IntStart, $R1IntStop, $R2IntStart, $R2IntStop);

}

sub getLastMapped {

	#return the position of the last mapped base relative to the read

	my ($cig, $ori) = @_;
	
	#reverse cigar if read is reverse
	if (($ori eq '-') or ($ori eq 'r')) { $cig = reverseCigar($cig); }
	
	
	#get all regions from cigar
	my (@letters, @numbers);
	getCigarParts($cig, \@letters, \@numbers);
	
	#remove any operations that don't consume query
	for my $i (0..$#letters) {
			if ($letters[$i] !~ /[MIS]/) {
			splice(@letters, $i, 1);
			splice(@numbers, $i, 1);
		}
	}
	
	my $lastMapped;
	for my $i (reverse(0..$#letters)) { #look through CIGAR from start to end
		if ($letters[$i] =~ /M/) { #if we've found the first match
			#lastMapped is the sum of the operations up to and including this one
			$lastMapped = eval join("+", @numbers[0..($i)]);		
			last;
		}
	}

	return $lastMapped;

}
