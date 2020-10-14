#!/usr/bin/perl
# It's not a bug, it's a feature

use strict;
use warnings;

use lib ".";

use ViralIntegration;
use Getopt::Long;

my $cutoff = 20; # default clipping cutoff
# mapped reads must have this number or more mapped bases
# unmapped reads must have fewer than this number of mapped bases
my $thresh = 0.95; #default amount of read that must be covered by alignments for rearrangement
my $tol = 5; #when processing CIGARS, combine any IDPN elements between matched regions with this number of bases or less
my $n_estimate_tlen = 10000; #  maximum number of pairs to use to estimate mean fragment length
my $viral;
my $host;
my $output = "integrationSites.txt";
my $bed;
my $merged;
my $verbose;

my $help;

our $warn_flag = 0;
$SIG{__WARN__} = sub { $warn_flag = 1; CORE::warn(@_) };

GetOptions('cutoff=i' => \$cutoff,
		   'viral=s'  => \$viral,
		   'host=s'  => \$host,
		   'output=s' => \$output,
		   'bed=s'    => \$bed,
		   'merged=s' => \$merged,
		   'tol=i'    => \$tol,
		   'n-estimate-tlen' => \$n_estimate_tlen,
		   'verbose'  => \$verbose,
		   'help',    => \$help);

if ($help) { printHelp(); }

unless ($viral and $host) { printHelp(); }

### hash arrays will store putative discordant pairs from both alignments
### make one array for R1 and R2, and one for viral and host (four total)
my %viralR1;
my %viralR2;
my %hostR1;
my %hostR2;

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

# sum template lengths from host and virus alignments, and use hash table to 
# keep track of reads that have already been added (to avoid double-counting)
my $virus_tlen;
my $host_tlen;
my %virus_tlen_reads;
my %host_tlen_reads;

# keep track of length of each reference for later when outputting locations in the host/virus genome
my %hostRlen;
my %virusRlen;

# to keep count of template lengths
my $counter = 0;

open (VIRAL, $viral) || die "Could not open viral alignment file: $viral\n";
if ($verbose) { print "Processing viral alignment...\n"; }
while (my $vl = <VIRAL>) {
	if ($vl =~ /^@/) { 
		if ($vl =~ /\@SQ/) {
		my @parts = split(' ', $vl);
		my ($ref) = ($parts[1] =~ /^SN:(.+)/);
		my ($len) = ($parts[2] =~ /^LN:(\d+)/);
		$virusRlen{$ref} = $len
		}
	next; 
	}

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
	if ($parts[1] & 0x100) { next(); } # skip alignments that are not primary
	
	# if read is proper pair, get template length
	if ($parts[1] & 0x2) {
		# if we haven't already counted enough reads
		if ($counter < $n_estimate_tlen) {
			unless (exists $virus_tlen_reads{$parts[0]}) {
				$virus_tlen += abs($parts[8]);
				$counter++;
				$virus_tlen_reads{$parts[0]} = '';
			}
		}
	}
	
	if ($parts[1] & 0x2) { next(); } # skip mapped in proper pair
	
	#for read pairs where one read is mapped and the other is not:
	#want reads where if not mapped, mate is mapped, or vice versa
	
	#but if allowing some clipping, might have case were one read has small portion mapped facing into the junction
	#so in this case also need to consider reads where both mates are mapped
	
	#so essentially only exclude reads where both reads are unmapped
	
	#if (($parts[1] & 0x8) and ($parts[1] & 0x4)) { next; }

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

# make sure that we found at least one reference length in the sam header, otherwise we don't have the information we'd need
unless (scalar keys %virusRlen > 0) {
	die "couldn't get length of reference viruses from sam header.  please check that header is intact";
}


# need to keep track of how many host alignments we used for estimating template length
my $virus_counter = $counter;
$counter = 0;

### Collect discordant reads from host genome
### This is the same process as for viral junctions
open (HOST, $host) || die "Could not open host alignment file: $host\n";
if ($verbose) { print "Processing host alignment...\n"; }
while (my $hl = <HOST>) {
	if ($hl =~ /^@/) { 
		if ($hl =~ /\@SQ/) {
		my @parts = split(' ', $hl);
		my ($ref) = ($parts[1] =~ /^SN:(.+)/);
		my ($len) = ($parts[2] =~ /^LN:(\d+)/);
		$hostRlen{$ref} = $len
		}
	next; 
	}

	my @parts = split(' ', $hl); #get fields from each alignment
	
	if ($parts[1] & 0x800) { next(); } # skip supplementary alignments
	if ($parts[1] & 0x100) { next(); } # skip alignments that are not primary
	
	# if read is proper pair, get template length
	if ($parts[1] & 0x2) {
		# if we haven't already counted enough reads
		if ($counter < $n_estimate_tlen) {
			unless (exists $host_tlen_reads{$parts[0]}) {
				$host_tlen += abs($parts[8]);
				$counter++;
				$host_tlen_reads{$parts[0]} = '';
			}
		}
	}
	
	if ($parts[1] & 0x2) { next(); } # skip mapped in proper pair
	
	unless ((exists $viralR1{$parts[0]}) and (exists $viralR2{$parts[0]})) { next; } #only consider reads from viral alignment
	
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
	if ($parts[1] & 0x40) 		{ $hostR1{$ID} = join("xxx", $seq, $seqOri, $ref, $start, $cig, $hSec, $hSup, ($editDist+$editDist2));}
	elsif ($parts[1] & 0x80)  	{ $hostR2{$ID} = join("xxx", $seq, $seqOri, $ref, $start, $cig, $hSec, $hSup, ($editDist+$editDist2)); }

}
close HOST;

# make sure that we found at least one reference length in the sam header, otherwise we don't have the information we'd need
unless (scalar keys %hostRlen > 0) {
	die "couldn't get length of host chromosomes from sam header.  please check that header is intact";
}

# need to keep track of how many host alignments we used for estimating template length
#first check that there was at least one proper pair of each
my $tlen;
my $host_pairs = keys %host_tlen_reads;
my $virus_pairs = keys %virus_tlen_reads;
# if we didn't find any proper pairs, then we can't estimate template length
# assume most pairs were merged, so set to 0
if ($host_pairs < 1 and $virus_pairs < 1) {
	$tlen = 0;
	print("no proper pairs found: setting template length to 0\n");
}
else {
	$tlen = ($host_tlen + $virus_tlen) / ($counter + $virus_counter);
	#https://stackoverflow.com/questions/178539/how-do-you-round-a-floating-point-number-in-perl
	$tlen = int($tlen + 0.5);
	print("template length estimated to be $tlen from host and virus alignments\n");
}


### Look for evidence of integration sites by identifying discordant read-pairs
### Need to compare the junction sites in the viral and host genomes to see if there is any overlap (ambiguous bases)
### Integration coordinates are not exact, but assigned to the middle part of the host read
my @outLines;
if ($verbose) { print "Finding discordant read-pairs...\n"; }
foreach my $key (keys %hostR1) {

	#check that read ID is found for all hashes
	unless ( exists $viralR1{$key} ) { next; } 
	unless ( exists $viralR2{$key} ) { next; }
	unless ( exists $hostR2{$key} ) { next; }

	# Find discordant read pairs
	my @intData = findDiscordant($key, $viralR1{$key}, $viralR2{$key}, $hostR1{$key}, $hostR2{$key}, $tlen, \%hostRlen, \%virusRlen);
	unless (@intData) { next; }	

	# Once the positions are collected, put together the output	
	my $outLine = join("\t", @intData);
	my $combReads = join("xxx", (split("xxx", $viralR1{$key}))[0], (split("xxx", $viralR2{$key}))[0]);
	if (defined $outLine) { push(@outLines, join("\t", ($outLine, $key, $combReads))); }
}



if ($verbose) { print "Writing output...\n"; }

my $header = "Chr\tIntStart\tIntStop\tVirusRef\tVirusStart\tVirusStop\tNoAmbiguousBases\tOverlapType\tOrientation\tHostSeq\tViralSeq\tAmbiguousSeq\tHostEditDist\tViralEditDist\tTotalEditDist\tPossibleHostTranslocation\tPossibleVectorRearrangement\tHostPossibleAmbiguous\tViralPossibleAmbiguous\tType\tReadID\tmerged\n";		
		
printOutput($output, $header, @outLines); #write to outfile: if no sites detected will be header only

if ($bed)    { printBed($bed, @outLines); }

if ($merged) { printMerged($bed, $merged); }

exit;


#----------------------------------------------------------------------------------------------------------------------------------------
### Subroutines ###
#----------------------------------------------------------------------------------------------------------------------------------------

sub findDiscordant {
### Check if there is overlap between the integration sites
### in host and viral integration sites
### returns integration start/stop relative to read sequence

### BWA Alignments are 1-based
	my ($key, $vR1, $vR2, $hR1, $hR2, $tlen, $hostRlen, $virusRlen) = @_;
	
	my ($seq1, $vR1ori, $vR1ref, $vR1start, $vR1cig, $vR1sec, $vR1sup, $vNM1) = (split('xxx', $vR1))[0, 1, 2, 3, 4, 5, 6, -1];
	my ($seq2, $vR2ori, $vR2ref, $vR2start, $vR2cig, $vR2sec, $vR2sup, $vNM2) = (split('xxx', $vR2))[0, 1, 2, 3, 4, 5, 6, -1];
	my ($hR1ori, $hR1ref, $hR1start, $hR1cig, $hR1sec, $hR1sup, $hNM1) = (split('xxx', $hR1))[1, 2, 3, 4, 5, 6, -1];
	my ($hR2ori, $hR2ref, $hR2start, $hR2cig, $hR2sec, $hR2sup, $hNM2) = (split('xxx', $hR2))[1, 2, 3, 4, 5, 6, -1];

	# discordant read-pair can be either one mate mapped and one unmapped
	# or one mate mostly mapped and the other with only a small mapped region
	
	# find if R1 and R2 are mostly mapped or unmapped
	my ($vR1map, $vR1mapBP) = isMapped($vR1cig, $cutoff);
	my ($vR2map, $vR2mapBP) = isMapped($vR2cig, $cutoff);
	my ($hR1map, $hR1mapBP) = isMapped($hR1cig, $cutoff);
	my ($hR2map, $hR2mapBP) = isMapped($hR2cig, $cutoff);
	
	# check that we do have complementary discordant read-pairs in host and viral alignments
	# check that both reads are not either both mapped or unmapped in each alignment
	unless (($vR1map ne $vR2map) and ($hR1map ne $hR2map)) { return; }
	unless (($vR1map eq $hR2map) and ($hR2map eq $vR1map)) { return; }
	
	# if a read-pair has gotten to this point, check that number of unmapped bases in 'map'
	# is less than $cutoff
	# do this to make sure we're reasonably confident of mapping and avoid reads with relatively
	# small mapped regions (eg 91S25M34S) getting through
	
	my $readlen1 = length($seq1);
	my $readlen2 = length($seq2);
	
	# get reference for virus and host based on which is mapped
	my ($vRef, $vSeq, $vCig, $vOri, $vPos, $vSec, $vSup, $vNM);
	my ($hRef, $hSeq, $hCig, $hOri, $hPos, $hSec, $hSup, $hNM);
	my ($intNM, $junct);
	
	#if host mapped R1, virus mapped R2
	if ($hR1map eq "map") { 
		#check unmapped bases is less than cutoff
		unless (($readlen1 - $hR1mapBP) < $cutoff) { return; } 
		unless (($readlen2 - $vR2mapBP) < $cutoff) { return; } 
		
		$hRef = $hR1ref; 
		$hSeq = $seq1; 
		$hNM = $hNM1;
		$hCig = $hR1cig;
		$hOri = $hR1ori;
		$hPos = $hR1start;
		$hSec = $hR1sec;
		$hSup = $hR1sup;
		
		$vRef = $vR2ref;
		$vSeq = $seq2;
		$vNM = $vNM2;
		$vCig = $vR2cig;
		$vOri = $vR2ori;
		$vPos = $vR2start;
		$vSec = $vR2sec;
		$vSup = $vR2sup;
		
		$intNM = $vNM + $hNM;
		
	}
	#if host mapped R2, virus mapped R1
	else { 
		
		#check unmapped bases is less than cutoff
		unless (($readlen2 - $hR2mapBP) < $cutoff) { return; } 
		unless (($readlen1 - $vR1mapBP) < $cutoff) { return; } 
		
		$hRef = $hR2ref; 
		$hSeq = $seq2; 
		$hNM = $hNM2;
		$hCig = $hR2cig;
		$hOri = $hR2ori;
		$hPos = $hR2start;
		$hSec = $hR2sec;
		$hSup = $hR2sup;
		
		$vRef = $vR1ref;
		$vSeq = $seq1;
		$vNM = $vNM1;
		$vCig = $vR1cig;
		$vOri = $vR1ori;
		$vPos = $vR1start;
		$vSec = $vR1sec;
		$vSup = $vR1sup;
		
		$intNM = $vNM + $hNM;
		
	}	


	# get properties of integration
	my ($hIntStart, $hIntStop, $vIntStart, $vIntStop, $isVecRearrange, $isHumRearrange, $isHumAmbig, $isVirAmbig, $nHAmbig, $nVAmbig, $hJunctSide, $vJunctSide);	
		
	#if host orientation is forward, junction is host/virus
	# otherwise it's virus/host
	if ($hOri eq 'f') {
		$junct = 'hv'; 
		$hJunctSide = 'left';
		# if the host and virus are mapped in the same orientation, need to reverse sides
		if ($vOri eq 'f') {
			$vJunctSide = 'left';
		}
		else {
			$vJunctSide = 'right';
		}
	}
	else {
		$junct = 'vh'; 
		$hJunctSide = 'right';
		if ($vOri eq 'f') {
			$vJunctSide = 'left';
		}
		else {
			$vJunctSide = 'right';
		}
	}

	# find if junction is host/virus or virus/host
	# find approximate location of junction on both host and virus side (assign to end of read closest to junction)
	# position indicated in SAM file is left-most mapping base: for forward reads this is 5' end, for reverse it's 3'
	
	# get length of reference to which read is aligned for host and virus
	my $hRefLen = $hostRlen->{$hRef};
	my $vRefLen = $virusRlen->{$vRef};
	
	($hIntStart, $hIntStop, $nHAmbig) = getIntPos($hCig, $hPos, length($hSeq), $vCig, length($vSeq), $hJunctSide, $tlen, $hRefLen);
		
	($vIntStart, $vIntStop, $nVAmbig) = getIntPos($vCig, $vPos, length($vSeq), $hCig, length($hSeq), $vJunctSide, $tlen, $vRefLen);

	# check for ambiguities
	if ((join(";", $vSup, $vSec)) eq "NA;NA") { $isVecRearrange = "no"; }
	else { $isVecRearrange = isRearrangeOrInt($vCig, $vOri, $vRef, $vPos, $vSup, $vSec, $vSeq, $thresh, $vNM, $intNM);}

	if ((join(";", $hSec, $hSup)) eq "NA;NA") { $isHumRearrange = "no"; }
	else { $isHumRearrange = isRearrangeOrInt($hCig, $hOri, $hRef, $hPos, $hSup, $hSec, $hSeq, $thresh, $hNM, $intNM);}
	 
	# check to see if location of host alignment is ambiguous: multiple equivalent alignments accounting for host part of read
	if ($hR1sec eq "NA") { $isHumAmbig = "no";}
	else { $isHumAmbig = isAmbigLoc($hOri, $hCig, $hSec, 'discordant', $hSeq, $tol);}
	
	# check to see if location of viral alignment is ambiguous: multiple equivalent alignments accounting for viral part of read
	if ($vSec eq "NA") { $isVirAmbig = "no";}
	else { $isVirAmbig = isAmbigLoc($vOri, $vCig, $vSec, 'discordant', $vSeq, $tol);}		
	


	return($hRef, $hIntStart, $hIntStop, $vRef, $vIntStart, $vIntStop, 
				"", 'discordant', $junct, $hSeq, $vSeq, '-', $hNM, $vNM, $intNM, 
				$isHumRearrange, $isVecRearrange, $isVirAmbig, $isHumAmbig, 'discordant');

}

sub printHelp {
	print "Pipeline for detection of viral integration sites within a genome\n\n";
	print "Usage:\n";
	print "\tperl discordant.pl --viral <sam> --host <sam> --cutoff <n> --output <out> --bed <bed> --tol <tol> --help\n\n";
	print "Arguments:\n";
	print "\t--viral:   Alignment of reads to viral genomes (sam)\n";
	print "\t--host:   Alignment of reads to host genome (sam)\n";
	print "\t--cutoff:  Minimum number of clipped reads to be considered (default = 20)\n";
	print "\t--tol:     Tolerance when combining short elements with neigbouring matched regions (default = 5)\n";
	print "\t--n-estimate-tlen: Number of pairs to use for estimating template length (default = 10000)\n";
	print "\t--output:  Output file for results (default = integrationSite.txt\n";
	print "\t--bed:     Print integrations sites to indicated bed file (default = NA)\n";
	print "\t--merged:  Merge bedfile into overlapping integration sites (default = NA)\n";
	print "\t--verbose: Print progress messages\n";
	print "\t--help:    Print this help message and quit\n";

	exit;
}

sub isMapped{

	my ($cig, $cutoff) = @_;
	
	# check if cigar is 'mapped' or not
	
	# read is considered to be unmapped if number of mapped bases is less than $cutoff
	# otherwise it's mapped ($cutoff or more mapped bases)
	
	# return 'map' or 'unmap' for mapped or unmapped read respectively
	# also return number of mapped bases
	
	if ($cig =~ /\*/) { return ('unmap', 0); }
	
	my (@letters, @numbers);
	getCigarParts($cig, \@letters, \@numbers);
	
	my $mappedBP = 0 ;
	for my $i (0..$#letters) {
		if ($letters[$i] =~ /M/) { $mappedBP += $numbers[$i]; }
	}

	if ($mappedBP >= $cutoff) { return ('map', $mappedBP); }
	return ('unmap', $mappedBP);
	
}

sub getIntPos {

	#my ($mapR1cig, $mapR1ori, $mapR1start, $readlen1, $mapR2cig, $mapR2ori, $mapR2start, $readlen2, $tlen, $junct) = @_;
	
	my ($cig, $pos, $rlen, $mateCig, $mateRlen, $junctionSide, $tlen, $reflen) = @_;
	
	#get positions of last mapped base facing in towards junction
	# output 0-based coordinates
	
	# give an estimate of integration site location based on template length/insert size
	
	# if the read is on the left side of the junction (ie the host read in a 'hv' juction
	# or the virus read in a 'vh' junction), then the start coordinate is the coordinate
	# of the base after the last mapped base in the read, 
	# and the stop coordinate is the hypothetical coordinate of the base 
	# before the first mapped base in the mate pair
	# this hypothetical coordinate is given by the sum of any clipped bases on the right
	# of the mapped right and the left of its mate, as well as the calculated insert length
	
	# for a read on the right side of the junction (ie the virus read in a 'hv' junction or
	# a host read in a 'vh' junction), then the start coordinate is the coordinate of the 
	# base after the hypothetical last mapped base of the mate, and the stop is the base before
	# the first mapped base of the read
	
	# cacluate insert length as template length minus sum of read lengths
	my $insert = $tlen - $rlen - $mateRlen;
	# insert length can't be less than 0
	if ($insert < 0) { $insert = 0; }
	
	my ($start, $stop, $ambig);
	
	if ($junctionSide eq 'left') {
		# get the number of mapped bases in the read
		my $rightClipped = $rlen - getLeftOrRightMapped($cig, 'right');
		
		my $mapped = $rlen - $rightClipped - (getLeftOrRightMapped($cig, 'left') - 1);
		
		# start is the last mapped base in the read
		$start = $pos + $mapped;
		
		# get number of left clipped bases in mate
		my $leftClippedMate = getLeftOrRightMapped($mateCig, 'left') - 1;
		
		# stop is start + insert + clipped bases
		$stop = $start + $rightClipped + $leftClippedMate + $insert;	
		
		# number of ambiguous bases is insert + clipped bases
		$ambig = $rightClipped + $leftClippedMate + $insert;
	}
	else {

		# stop is the base before the first mapped base
		$stop = $pos - 1;
		
		# get the number of soft-clipped bases on the left
		my $leftClipped = getLeftOrRightMapped($cig, 'left') - 1;
		
		# get the number of soft-clipped bases on the right of the mate
		my $rightClippedMate = $mateRlen - getLeftOrRightMapped($mateCig, 'right');
		
		# start is stop - insert - clipped bases
		$start = $stop - $insert - $leftClipped - $rightClippedMate;
		
		# number of ambiguous bases is insert + clipped bases
		$ambig = $insert + $leftClipped + $rightClippedMate;
		
	}
	
	# make sure that start isn't less than zero, and stop isn't more than the reference length
	if ($start < 0) { $start = 0; }
	if ($stop > $reflen) { $stop = $reflen; }

	return ($start, $stop, $ambig);
	

}

sub getLeftOrRightMapped {

	# return the position of the left-most or right-most mapped base (in the reference) 
	# relative to the read
	# output in 1-based numbering 
	# (ie in 150bp read, first base is 1 and last base is 150)
	
	# expected output:
	# getLeftOrRightMapped('1S3M2S', 'left') => 2
	# getLeftOrRightMapped('1S3M2S', 'right') => 4

	my ($cig, $type) = @_;
	
	# check that type is 'left' or 'right'
	if ($type ne 'left' and $type ne 'right') {
		die "type for getFirstOrLastMapped must be left or right"
	}
	
	#get all regions from cigar
	my (@letters, @numbers);
	getCigarParts($cig, \@letters, \@numbers);
	
	#remove any operations that don't consume query
	for my $i (reverse 0..$#letters) {
		if ($letters[$i] !~ /[MIS]/) {
			splice(@letters, $i, 1);
			splice(@numbers, $i, 1);
		}
	}
	
	my $mappedPos;
	
	# if we want the last mapped, we need to reverse letters and numbers backwards
	my @order;
	if ($type eq 'right') {
		@order = reverse(0..$#letters);
	}
	else {
		@order = (0..$#letters);
	}
	
	for my $i (@order) { #look through CIGAR from start to end
		if ($letters[$i] =~ /M/) { #if we've found the first match
			# last is the sum of the operations up to and including this one
			if ($type eq 'right') {
				$mappedPos = eval join("+", @numbers[0..($i)]);	
			}
			# first is the sum of operations up to but not including this one
			else {
				$mappedPos = eval join("+", @numbers[0..($i-1)]) || 0;
				$mappedPos++;
			}
			
			last;
		}
	}

	return $mappedPos;

}
