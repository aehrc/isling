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
my $tol = 5; #when processing CIGARS, combine any IDPN elements between M regions with this number of bases or less
my $viral;
my $human;
my $output = "short.txt";
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

	my ($cig, $combinedBases) = processCIGAR2($parts[5], $tol); # Process the CIGAR string to account for small insertions/deletions
	unless ($cig) { next; } # keep checking to make sure double clipped reads don't sneak through
	
	unless ($cig =~ /(^\d+[SH].*\d+[SH]$)/) { next; } # want alignments where both ends are clipped
	
	#number of mapped bases must be more than $cutoff
	my (@insert) = ($cig =~ /(\d+)M/g);
	my $mapBP = eval join("+", @insert);
	if ($mapBP < $cutoff) { next; }
	
	my $seq;
	my $dir;
	if ($parts[1] & 0x10) { 
		($seq) = reverseComp($parts[9]); 
		$dir   = 'r';
	}
	else				  { 
		$seq = $parts[9]; 
		$dir = 'f';
	}

	#get read ID, aligned viral reference, alignment start
	my ($ID, $ref, $start) = ($parts[0], $parts[2], $parts[3]);
	
	#get supplementary (SA) and secondary (XA) alignments in order to check for possible vector rearrangements
	my ($vSec, $vSup) = getSecSup($vl);
	my $editDist = getEditDist($vl);
	
	$viralIntegrations{join("xxx",($parts[0],$seq))} = join("xxx", $seq, $dir, $ref, $start, $cig, $vSec, $vSup, ($editDist+$combinedBases)); 
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

	my ($cig, $combinedBases) = processCIGAR2($parts[5], $tol); # Process the CIGAR string to account for complex alignments
	unless ($cig) { next; }
	
	unless ($cig =~ /(^\d+M.+\d+M$)/) { next; } #alignments must be matched at both ends
	unless ($cig =~ /I/) { next; } #must have inserted region
	
	#only want alignments where matched region is more than cutoff
	my ($match1, $match2) = ($cig =~ /^(\d+)M.+[ISDNPH](\d+)M$/);
	if (($match1 < $cutoff) or ($match2 < $cutoff)) { next; }
	
	#also want alignments where inserted region is more than cutoff
	#my (@insert) = ($cig =~ /(\d+)I/g);

	my $seq;
	my $Dir;
	if ($parts[1] & 0x10) { 
		($seq) = reverseComp($parts[9]); 
		$Dir   = 'r';
	}
	else				  { 
		$seq = $parts[9]; 
		$Dir = 'f';
	}
	
	#get read ID, aligned viral reference, alignment pos
	my ($ID, $ref, $start) = ($parts[0], $parts[2], $parts[3]);

	if (exists $viralIntegrations{join("xxx",($parts[0],$seq))}) { # only consider reads that were tagged from the viral alignment, no need to consider excess reads
	
		my ($hSec, $hSup) = getSecSup($hl);
		
		# get edit distance
		my $editDist = getEditDist($hl);
		
		# edit distance includes inserted bases, which aren't relevant for host alingment (these are the viral bases)
		# so remove the number of inserted bases from the edit distance
		my @letters;
		my @numbers;
		my $insertedBases = 0;
		getCigarParts($cig, \@letters, \@numbers);
		foreach my $i (0..$#letters) {
			if ($letters[$i] eq 'I') { $insertedBases += $numbers[$i]; }
		}
		
		$editDist -= $insertedBases;
		
		$humanIntegrations{join("xxx",($parts[0],$seq))} = join("xxx", $seq, $Dir, $ref, $start, $cig, $hSec, $hSup, ($editDist+$combinedBases));
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
		
		#get ID and seq from key
		my ($ID, $seq) = split('xxx', $key);
		
		#analyse short integration
		my ($int1Data, $int2Data) = analyseShort($viralIntegrations{$key}, $humanIntegrations{$key}, $ID, $seq);
		#push to array of data
		if ($int1Data) {
			push(@outLines, $int1Data);
			push(@outLines, $int2Data);
		}
	}
}


#print file
my $header = "Chr\tIntStart\tIntStop\tVirusRef\tVirusStart\tVirusStop\tNoAmbiguousBases\tOverlapType\tOrientation\tHostSeq\tViralSeq\tAmbiguousSeq\tHostEditDist\tViralEditDist\tTotalEditDist\tPossibleHostTranslocation\tPossibleVectorRearrangement\tHostPossibleAmbiguous\tViralPossibleAmbiguous\tType\tReadID\tmerged\n";		

if ($verbose) { print "Writing output...\n"; }

printOutput($output, $header, @outLines); #write to outfile: if no sites detected will be header only

if ($bed)    { printBed($bed, @outLines); }



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
	print "\t--thresh:	Amount of read that must be covered by one alignment to be considered rearrangement (default = 0.95)\n";
	print "\t--output:  Output file for results (default = short.txt\n";
	print "\t--bed:     Print integrations sites to indicated bed file (default = NA)\n";
	print "\t--merged:  Merge bedfile into overlapping integration sites (default = NA)\n";
	print "\t--verbose: Print progress messages\n";
	print "\t--help:    Print this help message and quit\n";

	exit;
}

sub analyseShort{

#find genomic coordinates of start and stop of each side of short integration

	my ($viralData, $hostData, $ID, $seq) = @_;
	
	#get data from input
	my ($vDir, $vRef, $vPos, $vCig, $vSec, $vSup, $vNM) = (split("xxx",$viralData))[1,2,3,4,5,6, -1];
	my ($hDir, $hRef, $hPos, $hCig, $hSec, $hSup, $hNM) = (split("xxx",$hostData))[1,2,3,4,5,6, -1];
	
	#arrays to store output
	#attributes of output:
	#0 - Host ref
	#1 - Host start
	#2 - Host stop
	#3 - Viral ref
	#4 - Viral start
	#5 - Viral stop
	#6 - no. ambigouous bases
	#7 - overlap type (gap, overlap, none)
	#8 - host seq (matched region in human alignment from end of read)
	#9 - viral seq (matched region in viral alignment)
	#10 - ambig seq (matched in both human and viral alignments)
	#11 - host secondary alignments
	#12 - viral secondary alignments
	#13 - possible translocation?
	#14 - possible vector rearrangement?
	#15 - host ambiguous location (is there a secondary alignment with same cigar as primary)
	#16 - viral ambiguous location (is there a secondary alignment with same cigar as primary)
	#17 - readID
	#18 - sequence
	
	#plus info about mate?  ie a way to link both integration sites from the one read
	#could look something like: a column with MATE: mateHStart;mateHStop;matevStart;matevStop;
	#or maybe don't need to include this, because can always pull out this information by filtering for both same ID and sequence
	
	#get int location relative to read
	
	#get start and stop locations of matched regions relative to read
	my @hMatched = getMatchedRegions($hCig, $hDir); #this should be start and end of read, plus some possible internal ones
	my @vMatched = getMatchedRegions($vCig, $vDir); #this should be middle of the read somewhere
	
	#get inserted region in viral read
	my $hInsert = getInserted($hCig, $hDir);
	
	#make sure that inserted bases in human alignment are in the same part of read as matched region in viral alignment
	#since matched regions can break up inserted region, if there is more than one inserted region, 
	#consider from the start of the first inserted region to the end of the last inserted region
	
	my $vMatchStart = (split('xxx', $vMatched[0]))[0];
	my $vMatchStop = (split('xxx', $vMatched[-1]))[1];
	
	my $hInsertStart = (split('xxx', $hInsert))[0];
	my $hInsertStop = (split('xxx', $hInsert))[1];
	
	#inserted region and matched region share bases if hStart < vStop and vStart < hStop
	unless (($vMatchStart < $hInsertStop) and ($hInsertStart < $vMatchStop)) { return; }
	
	#calculate overlap type for each
	#overlap is calculated based on matched regions
	
	my ($hMatch1Start, $hMatch1Stop) = split('xxx', $hMatched[0]);
	my ($hMatch2Start, $hMatch2Stop) = split('xxx', $hMatched[-1]);

	my $overlap1 = gapOrOverlap($hMatch1Stop, $vMatchStart);
	my $overlap2 = gapOrOverlap($vMatchStop, $hMatch2Start);
	
	#number of ambiguous bases
	#always calculate as right - left (relative to read)
	## ie for $ambig1, orientation is host-virus, so right is virus and left is host
	
	my $ambig1 = abs($hMatch1Stop - $vMatchStart) + 1; #1-based, so need to subtract 1
	my $ambig2 = abs($hMatch2Start - $vMatchStop) + 1;
	
	#number of inserted bases:
	my $inserted = $hInsertStop - $hInsertStart + 1;

	#get number of aligned bases for virus, first and second human alignment
	my ($hMatch1Len, $hMatch2Len, $vMatchLen);
	$hMatch1Len = $hMatch1Stop - $hMatch1Start + 1; #coords relative to read are 1-based
	$hMatch2Len = $hMatch2Stop - $hMatch2Start + 1; 
	$vMatchLen = $vMatchStop - $vMatchStart + 1; 
	
	#number of unambiguously aligned bases must be more than cutoff
	if ($overlap1 eq 'overlap') { unless (($hMatch1Len - $ambig1) > $cutoff) { return; } }
	if ($overlap2 eq 'overlap') { unless (($hMatch2Len - $ambig2) > $cutoff) { return; } }
	if (($overlap1 eq 'overlap') and ($overlap2 eq 'overlap')) {
		unless (($vMatchLen - $ambig1 - $ambig2) > $cutoff) { return; }
	}
	elsif ($overlap1 eq 'overlap') {
		unless (($vMatchLen - $ambig1) > $cutoff) { return; }
	}
	elsif ($overlap2 eq 'overlap') {
		unless (($vMatchLen - $ambig2) > $cutoff) { return; }
	}
	
	#calculate start and stop genomic positions  
	my ($hg1Start, $hg1Stop, $hg2Start, $hg2Stop, $vg1Start, $vg1Stop, $vg2Start, $vg2Stop);
	
	#get genomic coordinates of matched regions for host
	#note that if the read is reverse, start and stop will be swapped (so need to use start for first insertion and stop for second)
	#note that output of getGenomicCoords is 0-based
	my ($hgMatch1Start, $hgMatch1Stop) = getGenomicCoords($hMatch1Start, $hMatch1Stop, $hPos, $hDir, $hCig);
	my ($hgMatch2Start, $hgMatch2Stop) = getGenomicCoords($hMatch2Start, $hMatch2Stop, $hPos, $hDir, $hCig);
	
	#note that if the read is forward, (start, stop) will be ascending and opposite for reverse read
	#genomic coordinates refer to location of integration site (not mapped region)
	#in the case of an overlap, this is the genomic coordinates of the overlapped bases (which could be either host or viral)
	#in the case of a gap or a clean junction, this is the genomic coordinates of the last mapped base
	#don't give the coordinates of the gapped bases because they don't align to the host, and therefore their genomic coordinates are not defined
	
	##host coordinates
	if ($overlap1 eq 'overlap') { ($hg1Start, $hg1Stop) = ($hgMatch1Stop - 1, $hgMatch1Stop + $ambig1); }
	else 						{ ($hg1Start, $hg1Stop) = ($hgMatch1Stop - 1, $hgMatch1Stop); }
	#for second matched region, one or more deletions/matched regions in the middle of the read will mean that the coordinates 
	#aren't right next to first site
	#but still use second matched region
	if ($overlap2 eq 'overlap') { ($hg2Start, $hg2Stop) = ($hgMatch2Start - $ambig2, $hgMatch2Start + 1); }
	else 						{ ($hg2Start, $hg2Stop) = ($hgMatch2Start, $hgMatch2Start + 1); }		

	## viral coordinates
	my ($vgMatchStart, $vgMatchStop) = getGenomicCoords($vMatchStart, $vMatchStop, $vPos, $vDir, $vCig); #this assumes a single matched region, so could be problematic if there are more than one
	if ($overlap1 eq 'overlap') { ($vg1Start, $vg1Stop) = ($vgMatchStart - $ambig1, $vgMatchStart  + 1); }
	else 						{ ($vg1Start, $vg1Stop) = ($vgMatchStart, $vgMatchStart + 1); }
	if ($overlap2 eq 'overlap') { ($vg2Start, $vg2Stop) = ($vgMatchStop - 1, $vgMatchStop + $ambig2); }
	else 						{ ($vg2Start, $vg2Stop) = ($vgMatchStop - 1, $vgMatchStop); }
	
	#calculate total edit distance
	my $totalNM1 = $vNM + $hNM;
	my $totalNM2 = $vNM + $hNM;
	my $intNM = $vNM + $hNM;
	
	if ($overlap1 eq 'gap') { $totalNM1 += $ambig1; $intNM += $ambig1; }
	if ($overlap2 eq 'gap') { $totalNM2 += $ambig2; $intNM += $ambig2; }
	
	#check for ambiguity and rearrangement
	my $vRe;
	if ((join(";", $vSup, $vSec)) eq "NA;NA") 	{ $vRe = "no"; }
	else 										{ $vRe = isRearrangeOrInt($vCig, $vDir, $vRef, $vPos, $vSec, $vSup, $seq, $thresh, $vNM, $intNM);}

	my $hRe;
	if ((join(";", $hSup, $hSec)) eq "NA;NA") 	{ $hRe = "no"; }
	else 										{ $hRe = isRearrangeOrInt($hCig, $hDir, $hRef, $hPos, $hSec, $hSup, $seq, $thresh, $hNM, $intNM);}
	
	
	#check to see if location of human alignment is ambiguous: multiple equivalent alignments accounting for human part of read
	my $hAmbig;
	if ($hSec eq "NA") 	{ $hAmbig = "no";}
	else				{ $hAmbig = isAmbigLoc($hDir, $hCig, $hSec, 'short', $seq, $tol);}
	
	#check to see if location of viral alignment is ambiguous: multiple equivalent alignments accounting for viral part of read
	my $vAmbig;
	if ($vSec eq "NA") 	{ $vAmbig = "no";}
	else 				{ $vAmbig = isAmbigLoc($vDir, $vCig, $vSec, 'short', $seq, $tol);}

	#extract relevant parts of read
	## host sequence is only those bases that must come from the host - excluding any overlapped bases
	## similarly, viral sequence is only those bases that must come from the virus/vector
	## ambiguous bases 
	
	my ($hostSeq1, $ambigSeq1, $viralSeq, $ambigSeq2, $hostSeq2);	
	
	#first host sequence and ambiguous sequence
	if ($overlap1 eq 'overlap') {
		$hostSeq1 = substr($seq, $hMatch1Start - 1, $hMatch1Stop - $hMatch1Start + 1 - $ambig1);
		$ambigSeq1 = substr($seq, $hMatch1Stop - $ambig1, $ambig1);
	}
	else {
		$hostSeq1 = substr($seq, $hMatch1Start - 1, $hMatch1Stop - $hMatch1Start + 1);
		$ambigSeq1 = substr($seq, $hMatch1Stop, $ambig1);
	}
	
	#viral sequence
	if (($overlap1 eq 'overlap') & ($overlap2 eq 'overlap')) { #both sites have overlaps
		$viralSeq = substr($seq, $vMatchStart - 1 + $ambig1, $vMatchStop - $vMatchStart + 1 - $ambig1 - $ambig2 );
	}
	elsif ($overlap1 eq 'overlap') { #left site has overlap but not right
		$viralSeq = substr($seq, $vMatchStart - 1 + $ambig1, $vMatchStop - $vMatchStart + 1 - $ambig1);
	}
	elsif ($overlap2 eq 'overlap') { #right site has overlap but not left
		$viralSeq = substr($seq, $vMatchStart - 1, $vMatchStop - $vMatchStart + 1 - $ambig2);
	}
	else { #neither sites have overlaps
		$viralSeq = substr($seq, $vMatchStart - 1, $vMatchStop - $vMatchStart + 1);
	}
	
	#second host sequence and ambiguous sequence
	if ($overlap2 eq 'overlap') {
		$hostSeq2 = substr($seq, $hMatch2Start - 1 + $ambig2, $hMatch2Stop - $hMatch2Start + 1 - $ambig2);
		$ambigSeq2 = substr($seq, $hMatch2Start - 1, $ambig2);
	}
	else {
		$hostSeq2 = substr($seq, $hMatch2Start - 1, $hMatch2Stop - $hMatch2Start + 1);
		$ambigSeq2 = substr($seq, $hMatch2Start - 1 - $ambig2, $ambig2);
	}
	
	#output:
	## two arrays, one for each integration
	## hRef, hStart, hStop, vRef, vStart, vStop, noAmbigbases, overlapType, orientation, hostseq, viralseq, ambigSeq, hostSec, viralSec, possible translocation, possible vector rearrangement, host ambiguous, viral ambiguous, readID, seq
	
	my $int1Data = join("\t", $hRef, $hg1Start, $hg1Stop, $vRef, $vg1Start, $vg1Stop, $ambig1, $overlap1, 'hv', $hostSeq1, $viralSeq, $ambigSeq1, $hNM, $vNM, $totalNM1, $hRe, $vRe, $hAmbig, $vAmbig, 'short', $ID, $seq);
	my $int2Data = join("\t", $hRef, $hg2Start, $hg2Stop, $vRef, $vg2Start, $vg2Stop, $ambig2, $overlap2, 'vh', $hostSeq2, $viralSeq, $ambigSeq2, $hNM, $vNM, $totalNM2, $hRe, $vRe, $hAmbig, $vAmbig, 'short', $ID, $seq);
	
	return ($int1Data, $int2Data);
}

sub getInserted {

	#get inserted portion of read, based on CIGAR, and noting that some inserted regions may be broken up by matched regions
	#so if there is more than one inserted region, assume that inserted region stretches from start of first one to end of last one
	
	my ($cig, $dir) = @_;
	
	if (($dir eq '-') or ($dir eq 'r')) { $cig = reverseCigar($cig); }
	
	#get all regions from cigar
	my (@letters, @numbers);
	getCigarParts($cig, \@letters, \@numbers);
	
	#in order to get correct coordinates relative to read (query)
	#need to remove any cigar operations that don't consume query (D, N, H, P)
	for my $i (0..$#letters) {
		if ($letters[$i] =~ /[DNHP]/) {
			splice(@letters, $i, 1);
			splice(@numbers, $i, 1);
		}
	}
	
	#for each match in @letters, get position in read where inserted region starts and ends
	my ($start, $matchedLen, $end, @aligns);
	for my $i (0..$#letters) { 
		if ($letters[$i] eq "I") { 
			#get matched length
			$matchedLen = $numbers[$i];
		
			#get start
			$start = 1;
			for my $j (0..($i-1)) { $start += $numbers[$j]; } #sum up numbers up to where match starts
				
			#get end
			$end = ($start + $matchedLen - 1);
				
			#append to @aligns
			push(@aligns, "${start}xxx${end}");
		
			}
		}

	#get start of first insert
	my $iStart = (split('xxx', $aligns[0]))[0];
	
	#get end of last insert
	my $iStop = (split('xxx', $aligns[-1]))[1];
	
	return (join('xxx', $iStart, $iStop));

}
