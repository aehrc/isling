#!/usr/bin/perl

use strict;
use warnings;

package ViralIntegration;
use Exporter;
our @ISA = ('Exporter');
our @EXPORT = qw(isRearrangeOrInt getEditDist processCIGAR processCIGAR2 extractSeqCoords extractCoords extractSequence gapOrOverlap getCigarParts getGenomicCoords getMatchedRegion getMatchedRegions getSecSup isAmbigLoc isRearrange printOutput printBed printMerged reverseComp reverseCigar);


##### subroutines #####

sub extractSeqCoords {
### Extract coordinates of the human/viral sequences relative to read
### these are the bases we know must come from human/virus (so not including any overlapped bases)

	my ($ori, $align, $overlap, $length, $overlaptype) = @_;

	if ($ori eq "+") {  #mapped part is first
		if ($overlaptype eq "overlap") { return (1, $align - $overlap); }
		else { return (1, $align); }
	} 					
	else {  #mapped part is second
		if ($overlaptype eq "overlap") { return (($length - $align + $overlap + 1) , $length); }
		else { return (($length - $align + 1), $length); }
	}
}

sub extractCoords {
### Extract 0-based genomic coordinates for matched region based on start and stop relative to read
	my ($alig, $overlap, $overlaptype, $pos, $ori) = @_;
	
	#$alig is number of aligned bases 
	#$overlap is number of ambigous bases (gap or overlap)
	#$overlaptype is type of overlap (gap, overlap, none)
	#$pos is 1-based mapping coordiate of first mapped base (from BWA)
	#$dir is 'f' for reads mapped in the forward direction and 'r' for reads mapped in the reverse direction
	my ($intGStart, $intGStop);
	
	if ($ori eq '+') {
		$intGStart = $pos + $alig - 1;
		$intGStop = $pos + $alig;
		if ($overlaptype eq "overlap") { $intGStart -= $overlap; }
	}
	else {
		$intGStart = $pos - 1;
		$intGStop = $pos;
		if ($overlaptype eq "overlap") { $intGStop += $overlap; }
	}
	
	return ($intGStart, $intGStop);
	
}

sub gapOrOverlap {
	
	#calculate if there's a gap or overlap based on the start of one matched region and the end of the previous region
	# pass in stop of first region and start of second region relative to read	
	
	my ($stop1, $start2) = @_;
	
	if ($stop1 + 1 == $start2) { return "none"; } #1-based, so if they're bookended then there's no gap
	if ($stop1 + 1 > $start2) { return "overlap"; }
	if ($stop1 + 1 < $start2) { return "gap"; }
	
	return;

}

sub getCigarParts {
	#split up a cigar string into two different arrays containing the letters and the numbers
	#return references to the two arrays
	
	my ($cigar, $let, $num) = @_;
	
	#get and reverse letters and numbers from cigar
	@$let = split(/\d+/, $cigar);
	@$num = split(/[A-Z]/, $cigar);

	#since numbers precede letters, always get one extra empty element in letter array (which is at start)
	shift @$let;
	
}

sub getGenomicCoords {
	#get genomic coordinates (0-based) of cigar operation (specified by start and stop relative to the read)
	#using cigar, direction and POS (1-based leftmost mapping position of first CIGAR operation that consumes a reference base [M, D, N])
	#need to consider CIGAR operations that consume read [MIS=X] or not [DNHP], 
	#and those that consume reference [MDN] or not [ISHP]
	
	#without processing the cigar, might have more than one matched region per alignment so want to get coordinates for one that is specified
	
	my ($start, $stop, $pos, $sense, $cig) = @_;
	
	#need to know what the coordinate of the start of the read is, given the left-most mapping position of M, D or N cigar operation
	
	#get all regions from cigar
	my (@letters, @numbers);
	getCigarParts($cig, \@letters, \@numbers);
	
	#need two new versions of the @numbers array:
	#version @rNumbers retains numbers for all @letters that consume query [MIS]
	#version @gNumbers retains numbers of all @letters that consume reference [MDN]
	#need these to make calculation of read- and genome-relative positions for each CIGAR operation
	my @rNumbers = @numbers;
	my @gNumbers = @numbers;
	
	for my $i (0..$#numbers) {
		#if doesn't consume query
		if ($letters[$i] !~ /[MIS]/) {
			$rNumbers[$i] = 0;
		}
		#if doesn't consume reference
		if ($letters[$i] !~ /[MDN]/) {
			$gNumbers[$i] = 0;
		}
	}
	
	
	#make arrays with info about start and stop position relative to read of each CIGAR operation
	#also genomic position (start and stop) for each CIGAR operation
	my ($rStart, $rStop, $gStart, $gStop);
	
	#loop over elements of CIGAR to calculate start and stop positions relative to read
	for my $i (0..$#numbers) {
		#add to rStart and rStop 
		#rStart and rStop are 1-based
		#if forward, need to sum from start of array to $i
		if (($sense eq 'f') or ($sense eq '+')) {
			#rStart for this position is 1 + (sum of @numbers up to but not including this position)
			$rStart = 1 + eval join("+", @rNumbers[0..($i-1)]);
			#rStop for this position is (sum of @numbers up to and including this position)
			$rStop = eval join("+", @rNumbers[0..$i]);
		}
		#read is reverse, need to sum from $i to end of array
		else {
			#rStart is 1+ sum from next element to end of array
			$rStart = 1 + eval join("+", @rNumbers[$i+1..$#rNumbers]);
			#rStop is sum from this element to end of array
			$rStop = eval join("+", @rNumbers[$i..$#rNumbers]);
			#gStart is 
		}
		#if we've found the right operation
		if (($start == $rStart) and ($stop == $rStop)) {
			#calculate gStart and gStop
			#coordinates relative to read are 1-based, but output 0-based 
			
			if ($sense eq 'f') {
				$gStart = $pos - 1 + eval join("+", @gNumbers[0..($i-1)]);
				$gStop = $pos + eval join("+", @gNumbers[0..$i]);
			}
			else {
				$gStart = $pos + eval join("+", @gNumbers[0..$i]);
				$gStop = $pos - 1 + eval join("+", @gNumbers[0..($i-1)]);
			}
			return  ($gStart, $gStop) ;
			
		}
		
	}
}

sub getMatchedRegion {
	#get the matched region from a cigar string, assuming it's at the start or end of the read
	
	my ($cigar, $dir) = @_;

	if (($dir eq '-') or ($dir eq 'r')) { $cigar = reverseCigar($cigar); }
	
	my ($end, $align);
	
	if	  (($cigar =~ /^(\d+)M/)) { $end = 'start';  ($align) = ($cigar =~ /^(\d+)M/);}
	elsif (($cigar =~ /(\d+)M$/)) { $end = 'end'  ;  ($align) = ($cigar =~ /[A-Z](\d+)M$/);}
	
	return ($end, $align);

}

sub getMatchedRegions {
	#get all the matched regions from a cigar string
	#don't assume that cigar has been simplified in any way
	#returns a 1D array of matched regions in the form startxxxend	(relative to read)

	my ($cigar, $dir) = @_;
	
	
	if (($dir eq '-') or ($dir eq 'r')) { $cigar = reverseCigar($cigar); }
	
	#check for no alignments
	unless ($cigar =~ /M/) { return; }
	
	#if whole cigar is matched
	if ($cigar =~ /^(\d+)M$/) { my ($match) = ($cigar =~ /^(\d+)M$/); return "1xxx${match}"; }
	
	#get letters and numbers from cigar
	my (@letters, @numbers);
	getCigarParts($cigar, \@letters, \@numbers);
	
	
	#in order to get correct coordinates relative to read (query)
	#need to remove any cigar operations that don't consume query (D, N, H, P)
	
	for my $i (reverse(0..$#letters)) {
		if ($letters[$i] =~ /[DNHP]/) {
			splice(@letters, $i, 1);
			splice(@numbers, $i, 1);
		}
	}

	#for each match in @letters, get position in read where alignment starts and length of alignment
	my ($start, $matchedLen, $end, @aligns);
	for my $i (0..$#letters) { 
		if ($letters[$i] eq "M") { 
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
	return @aligns;
}

sub getSecSup {
#get secondary (XA) and supplementary (SA) alignments from line

#note that secondary alignemnts have the form /XA:Z:(chr,pos,CIGAR,NM;)*/
#supplementary alignments have the form /SA:Z:(chr,pos,strand,CIGAR,PHRED,NM;)*/

	my ($line) = @_;
	
	#get supplementary alignments from SA field
	my ($sup, $sec);
	if ($line =~ /SA:.:.*;/) 	{ 
		($sup) = ($line =~ /SA:.:(.*?);\s/);
	
	}
	else 						{ $sup = "NA"; }	
	
	#Add get secondary alignment information from XA field
	if ($line =~ /XA:.:.*;/) 	{ 
		($sec) = ($line =~ /XA:.:(.*);\s/); 
		
		#two different fields have different formats: 
		#SA has /(chr,pos,strand,CIGAR,PHRED,NM;)*/ but XA has only /(chr,pos,CIGAR,NM;)*/ where pos starts with + or - to indicate strand
		#convert all XA alignments to SA form to make them consistent
		my $newSecs = "";
		my @secs = split(";", $sec); #split alignments into array
		my $secAlignment;
		foreach $secAlignment (@secs) {
			my ($chr, $pos, $cig, $nm) = (split(",", $secAlignment))[0,1,2,3]; #split alignment into components
			my ($strand, $newPos) = ($pos =~ /([\+\-])(\d+)/); #split strand and position info from orignal $pos
	
			#join elements to make new alignment string
			#note that XA doesen't contain phred info, so just put 0 here
			my $newSec = join(",", $chr, $newPos, $strand, $cig, "0", $nm);
			
			$newSecs = join(";", $newSec, $newSecs)
		}
		$sec = $newSecs;

	}
	else 	{ $sec = "NA"; }
	
	return ($sec, $sup);

}

sub getEditDist {
#get edit distance (NM)  from line

	my ($line) = @_;
	
	#get supplementary alignments from SA field
	my ($nm);
	($nm) = ($line =~ /NM:i:(\d+)\s/);

	return ($nm);
	
}

sub isAmbigLoc {
### check for ambiguity in location of integration site in human or viral genome
## look to see if there are any secondary alignments that are equivalent to primary alignment
## there is ambiguity if the CIGAR string is identical for the primary and a secondary alignment

	#$dir is direction of primary alignment ('f' or 'r')
	#$cig is cigar of primary alignment
	#$sec is secondary alignments in the form /(chr,pos,CIGAR,NM;)*/
	#$type is type of integration - either 'soft', 'discordant' or 'short' (need this to know how to process CIGAR)
	#$seq is the read sequence (need for processCIGAR for soft-clipped reads)
	#tol is max number of bases to combine for processCIGAR2
	my ($dir, $cig, $sec, $type, $seq, $tol) = @_;
	
	#check type is one of the three options
	unless (($type eq 'short') or ($type eq 'discordant') or ($type eq 'soft')) {
		print "wrong type input when checking for ambiguous location\n";
		return;
	}
	
	
	#if alignment is in reverse orientation, convert to forward by reversing the cigar
	if ($dir eq 'r') { $cig = reverseCigar($cig); }
	
	my $isAmbiguous = "";
	my $dummy;
	my @secs = split(";", $sec); # split secondary alignments to check one by one
	foreach my $secAlignment (@secs) {
	
		if ($secAlignment eq "NA") { next; } #check for no secondary alignments
		my ($secSense, $secCigar) = (split(",",$secAlignment))[2,3]; #get sense and cigar from 
		
		#process cigar as appropriate for read type
		my $secCigar1;
		if ($type eq 'soft') { ($secCigar1, $dummy) = processCIGAR($secCigar,$seq);}
		($secCigar, $dummy) = processCIGAR2($secCigar, $tol);
	
		if ($secSense eq '-') { $secCigar = reverseCigar($secCigar); } #reverse cigar if necessary
		
		#if cigars match, there is ambiguity
		if ($secCigar eq $cig) { $isAmbiguous = "yes"; } 
		if ($secCigar1) {
			if ($secCigar1 eq $cig) { $isAmbiguous = "yes"; }
		}
	}
	#if haven't found any ambiguity yet
	if ($isAmbiguous eq "") { $isAmbiguous = "no"; }
	
	return $isAmbiguous;

}

sub isRearrange {

	#check for possible vector rearrangements:
	#get start and end positions of all alignments: primary, supplementary, secondary
	#sort alignments and check for gaps between start of read, each alignment, and end of read
	#check if number of bases in gaps is less than $thresh * readlength
	
	my ($pCig, $pDir, $pRef, $pPos, $sup, $seq, $thresh, $pNM) = @_;
	# $pCig is primary alignment of the first reference
	# $pDir is the direction of the primary alignment 
	# $sup are the secondary and supplementary alignments
	# $thresh is the fraction of the read that must be accounted for by alignments to be considered a rearrangement
	
	
	
	#get location in read of matched regions and their lengths
	#make array of start of matched regions and their lengths in format startxxxend
	#need to take into account direction of read - convert everything to forward orientation
	#need to zeropad start and length so that each have the same number of digits, for sorting later
	my $readlen = length($seq);

	
	#note - currently considering all alignments together regardless of their reference: consider changing to check only that alignments to a sigle reference account for the whole read
	#for vector/virus we would want to do this because probably don't have different viruses recombining?
	#for human/host we would not want to do this could have translocations
	#could do this by making one @aligns per reference, and store all 
	#note - using 1-based indexing of read alignments
	
	my @aligns; #to store start and end of each alignment
	
	#first do primary alignment
	my @pAligns = getMatchedRegions($pCig, $pDir);
	foreach my $align (@pAligns) { push(@aligns, join("xxx", $align, $pRef, $pPos, $pDir, $pCig, $pNM)); }
	
	#then do rest of the alignments
	my @supAligns = split(";", $sup);
	foreach my $supAlign (@supAligns) {
	
		if (($supAlign eq "NA") or ($supAlign eq "")) { next; } #check for no secondary alignments
		my ($supRef, $supPos, $supSense, $supCig, $supNM) = (split(",",$supAlign))[0,1,2,3,-1]; #get info about this alignment
		
		#need to keep the other alignment info for later output
		my @curAligns = getMatchedRegions($supCig, $supSense);
		foreach my $align (@curAligns) { push(@aligns, join("xxx", $align, $supRef, $supPos, $supSense, $supCig, $supNM)); }
		
	}
	
	#sort array
	my @sorted = sort(zeroPad(@aligns));
	
	#remove any nested alignments
	#these are alignments that are completely encompassed by another alignment 
	#ie 001xxx100 and 010xxx020 are nested and 010xxx020 should be removed
	#also 001xxx009 and 001xxx010 are nested and 001xxx009 should be removed
	#also remove any equivalent alignments - these are alignments with identical start and stop postions
	
	#strategy: start at end of array and check if the current alignment and the previous alignment are nested
	#if so, remove the nested one (which will always be the current alignment, since the list is sorted)

	my $lastKept = $#sorted;
	for my $i (reverse((1..$#sorted))) {
		#get current and previous start and end
		my ($cS, $cE) = (split('xxx', $sorted[$i]))[0,1];
		my ($pS, $pE) = (split('xxx', $sorted[$i-1]))[0,1];
		
		#if two (or more) alignments are equivalent, remove the one with the smaller edit distance
		if (($cS == $pS) and ($cE == $pE)) {
			#get edit distances
			my $cNM = (split('xxx', $sorted[$i]))[-1];
			my $pNM = (split('xxx', $sorted[$i-1]))[-1];
			if ($pNM > $cNM) { #if previous is more than current, keep current and delete previous
				#first move current alignment to the position of the previous alignment, then delete current
				$sorted[$i-1] = $sorted[$i];
				splice(@sorted, $i, 1);
			}
			else { #otherwise, keep previous and delete current
				splice(@sorted, $i, 1); 	
			}
			next;
		}
		#if current alignment and previous alignment are nested, and the current alignment is the nested alignment
		if ((($cS > $pS) and ($cE < $pE)) or (($cS > $pS) and ($cE == $pE))) { 
			#just remove the one alignment (the current alignment)
			splice(@sorted, $i, 1); 
			next;
		}
		#if current alignment and previous alignment are nested, and the previous alignment is the nested alignment
		#eg current = 001xxx120, previous = 001xxx100
		if (($cS == $pS) and ($cE > $pE)) {
			#first move current alignment to the position of the previous alignment, then delete
			$sorted[$i-1] = $sorted[$i];
			splice(@sorted, $i, 1);
		}

		
	}	
	
	#check if alignments can account for whole read
	#check all alignments for gaps
	#to store the number of total bases gapped and number of gaps
	my ($gapBP, $gaps);
	$gapBP = (split('xxx', $sorted[0]))[0] - 1; #initialise to gap between start of read and start of first alignment
	if ($gapBP == 0) { $gaps = 0; }
	else { $gaps = 1; }
	
	#check all alingments for gaps
	for my $i (1..$#sorted) {
		#if at end of the alignments, check gap between end of read and last alignment
		if ($i == $#sorted) {
			#check if gap between end of read and last alignment
			unless ((split('xxx', $sorted[$i]))[1] == $readlen) {
				#increment gapBP and gap
				$gapBP += $readlen - (split('xxx', $sorted[$i]))[1];
				$gaps += 1;
			}
		} 
		#if start of current alignment is greater than or equal to (end of previous alignment + 1)
		unless ((split('xxx', $sorted[$i]))[0] <= (split('xxx', $sorted[$i-1]))[1] + 1 ) {
			#add number of gap bases between start of current alignment and end of previous alignment to total gapped bases
			 $gapBP +=  (split('xxx', $sorted[$i]))[0] - ((split('xxx', $sorted[$i-1]))[1] + 1);
			$gaps += 1
		
		}
	}
	
	#get total edit distance from all alignments
	my $NM = $gapBP;
	for my $align (@sorted) {
		$NM += (split("xxx", $align))[-1];
	}

	
	my $isRearrange; #to store result, "yes" or "no"
	
	if ($gapBP < ($readlen - ($thresh * $readlen)) and ($#sorted > 0)) { $isRearrange = "yes"; }
	else { $isRearrange = "no"; }
	
	
	#get start and stop 
	
	return ($isRearrange, $gapBP, $gaps, $NM, @aligns);
	
	
}

sub isRearrangeOrInt {
	#decide if a given read is more likely to result from an integration or a rearrangement
	#do this by comparing a 'total edit distance'
	
	#for integrations, the total edit distance is the sum of the virus and host alignment edit distances, 
	#plus the number of bases in the gap if there is one
	
	#for rearrangements, the total edit distance is the sum of all the edit distances of the alignments used to account for the bases in the read
	#plus the number of bases unaccounted for 
	
	#return 'yes' if the read is more likely a rearrangement
	#or 'no' if more likely an integration
	
	my ($cig, $dir, $ref, $pos, $sec, $sup, $seq, $thresh, $pNM, $intNM) = @_;
	
	my ($rearrangeNM, $isRearrange);
	
	#if there's no sec or sups just make rearrangeNM the total length of the read
	if ((join(";", $sup, $sec)) eq "NA;NA") { $rearrangeNM = length($seq); }
	
	#otherwise get it from isRearrange
	else { $rearrangeNM = (isRearrange($cig, $dir, $ref, $pos, (join(";", $sup, $sec)), $seq, $thresh, $pNM))[3];}
	
	#decide if read is more likely a rearrangement or integration based on edit distance
	#if edit distances are the same, err on the side of caution and assign it as a rearrangement
	if ($intNM < $rearrangeNM) 	{ $isRearrange = "no"; }
	else						{ $isRearrange = "yes"; }
	
	return $isRearrange;
	
}

sub printOutput {
### Print output file
	my ($outFile, $header, @outLines) = @_;	
	open (OUTFILE, ">$outFile") || die "Could not open output file: $outFile\n";
	print OUTFILE $header;
	foreach my $line (@outLines) { print OUTFILE "$line\n"; }
	close OUTFILE;
}

sub printBed {
### Print Bed output
### Bed starts are 0-based and their ends are 1-based
### 1	0	1 <- referes to the first base on chromosome 1

	my ($bedFile, @outLines) = @_;
	open (BEDFILE, ">$bedFile") || die "Could not open bedfile: $bedFile\n";
	foreach my $line (@outLines) {
		my @parts = split("\t",$line);
		print BEDFILE "$parts[0]\t$parts[1]\t$parts[2]\t$parts[-2]\t.\t$parts[8]\n";
	}
	close BEDFILE;
}

sub printMerged {
### merge bed into overlapping regions
	my ($bed, $merged) = @_;
	my $sorted = `sort -k1,1 -k2,2n $bed`;

	my $uni = int(rand(1000000));
	
	my $temp = join("", ("temp.", $uni, ".bed"));
	open (TEMP, ">$temp");
	print TEMP $sorted;
	close TEMP;

	my $mergeOut = `bedtools merge -c 5,4 -o count,distinct -i $temp`;
	open (MERGED, ">$merged") || die "Could not open merged file: $merged\n";
	print MERGED $mergeOut;
	close MERGED;
	
	unlink($temp);
}

sub processCIGAR {
	my ($oriCig,$seq) = @_;

	my ($clipped, $other, $order);

	# Order:
	#	1 - clipped,matched
	#	2 - matched,clipped

	if ($oriCig =~ /^\d+[SH].+\d+[SH]$/) { return; } # skip if a double clipped read got through
	elsif ($oriCig =~ /^(\d+)[SH](.+)$/) { ($clipped, $other, $order) = ($1,$2,1); }
	elsif ($oriCig =~ /^(.+[MIDP])(\d+)[SH]$/) { ($other, $clipped, $order) = ($1,$2,2); }
	# Note I don't save if it's hard/soft clipped, only the number of bases

	# The CIGAR String can have multiple combinations of the following features:
	# M - matched
	# S - soft clipped (bases kept in the read)
	# H - hard clipped (bases removed from the read)
	# I - inserted seqeunce
	# D - deleted sequence
	# P - padding

	# I, D and P can break up a matched segment (e.g. MIPM) but should still be counted as a single matched sequence
	# I needs to be included in the final "matched" count - because everything is done relative to the read
	# P needs to be removed in the final "matched" count  - because it is just adding space to help the read align
	# D needs to be removed in the final "matched" count  - because it adds space where base was missing

	# The full "matched" sequence, can be calculated as length(read) - length(clipped).
	# This will give the total amount of the read that was aligned (matched) and the part of the read that wasn't (clipped)
	# The assumption is that everything else is either matched or inserted (and therefore should still be counted)

	my $matched = length($seq) - $clipped;	
	
	# also calculate the number of bases added to the matched region, to add to the edit distance
	#get the number of bases originally in the matched region
	my @matchedRegions = getMatchedRegions($oriCig, '+');
	my $sum = 0;
	foreach my $region (@matchedRegions) {
		my ($start, $end) = split("xxx", $region);
		$sum += $end - $start + 1;
	}
	# number of bases added is final matched bases - original matched bases
	my $addedBases = $matched - $sum;
	
	# Construct simplified CIGAR string (only contains matched and clipped)
	if ($order == 1) { return(join('',($clipped,"S",$matched,"M")), $addedBases); }
	else		 { return(join('',($matched,"M",$clipped,"S")), $addedBases); }
}

sub processCIGAR2 {

	# The CIGAR String can have multiple combinations of the following features:
	# M - matched
	# S - soft clipped (bases kept in the read)
	# H - hard clipped (bases removed from the read)
	# I - inserted seqeunce
	# D - deleted sequence
	# P - padding
	# N - skipped region from reference (only relevant for introns)
	
	#process CIGAR to remove small I, D, P and N cigar operations that break up a matched region
	#total number of bases in small I, D, P and N operations should be less than $tol in order to be removed
	
	#don't assume that result of processing should be read with only one soft-clipped and one matched region 
	#(as was assumed by original processCIGAR subroutine)
	
	#CIGAR operations that consume query bases need to be incorporated into the nearest matched region ()
	#CIGAR operations that don't consume query bases can be discarded
	
	#strategy for processing:
	#if CIGAR contains no I, D or P, N operations just return orignal CIGAR
	#similarly, if there's no matched regions, return the original CIGAR
	
	#if there's only one matched region in the read, return the original CIGAR
	
	#otherwise consider each CIGAR operation - split into letters and numbers
	#if the CIGAR contains one or more [IDPN] operation between two M operations,
	#and the number of bases in those [IDPN] operations is less than $tol,
	#combine those [IDPN] operations with the two adjacent M operations to leave only 
	#one matched region
	
	#I operation consumes query, so these bases should be added to the number of bases in the two adjacent M operations
	#[DNP] operations don't consume query, so these bases can be discarded
	
	#also consider [IDPN] between a matched region and the start or end of the read
		#combine with the matched region
	#also consider [IDPN] between a matched region and a soft-clipped region (this shouldn't happen???)
	
	my ($oriCig, $tol) = @_;

	#if no [IDP] operations, return original CIGAR
	#or if doesn't contain any match operations, also return original CIGAR 
	#(this could happen if read is unmapped, so CIGAR is *)
	if ($oriCig !~ /[IDPNM]/) { return $oriCig, 0; }
	
	#get letters and numbers in CIGAR as array
	my (@letters, @numbers);
	getCigarParts($oriCig, \@letters, \@numbers);
	
	#if CIGAR only has one operation, return original CIGAR
	if ($#letters == 0) { return $oriCig, 0; }
	
	#add zero-length matched regions at start and end of letters and numbers in order to 
	#properly combine [IDPN] operations that occur at the start or end of read
	unshift(@letters, 'M');
	unshift(@numbers, '0');
	push(@letters, 'M');
	push(@numbers, '0');
	
	#look backwards through @letters, and merge any [IDPN] operations with a length of less than $tol
	#that are sandwiched between two M operations
	my $lastM = $#letters;
	my $bases;
	my $combinedBases = 0;
	for my $i (reverse((0..$#letters))) {
		#check if this element is M
		if ($letters[$i] =~ /M/) {
		
			#don't want to combine soft-clips - it is common to have a few bases soft-clipped at either end of the read
			#the ends of the read is the only place this would occur - from SAM format:
			#"S may only have H operations between them and the ends of the CIGAR string"
			
			#so if there's an S in the operations sandwiched between two M's
			#just move on
			my @Sops = grep { $_ =~ /S/ } @letters[$i..$lastM];
			if( $#Sops >= 0) { $lastM = $i; next; }
			
			#check the number of bases between this M and the last M
			
			#if there are no elements between $i+1 and $lastM-1, the number of bases is 0
			$bases = eval join "+", ( @numbers[$i+1..$lastM-1] || (0) ) ;
			if (($bases <= $tol) && ($lastM != $i)) { 
					#combine the elements between the two Ms
					#need to consider whether to add bases to total or not
					#only add bases to total if CIGAR operation consumes query ([MI])
					#otherwise discard bases
					my $total;
					for my $j ($i..$lastM) {
						if ($letters[$j] =~ /[MIS]/) { $total += $numbers[$j]; }
					}
				#add number of bases to count and cut out of cigar
				$combinedBases += $numbers[$i+1..$lastM-1] || 0;
				splice(@letters, $i+1, $lastM-$i);
				splice(@numbers, $i, $lastM-$i+1, $total);
			}
			#reassign $lastM for next step
			$lastM = $i;
		}
	}	
	
	#before reconstructing CIGAR, need to remove any 0M operations added
	for my $i (reverse((0..$#letters))) {
		if ($numbers[$i] == 0) {
			splice(@letters, $i, 1);
			splice(@numbers, $i, 1);
		}	
	}
	
	#reconstruct CIGAR
	my $newCig;
	foreach my $number (@numbers) {
		$newCig .= ($number.shift(@letters));	
	}
	
	#return both new cigar and number of bases converted to matched (for edit distance)
	
	return $newCig, $combinedBases;
}

sub reverseComp {
### return reverse compliment of sequence
	my ($seq) = @_;
	my $rev = reverse $seq;
	$rev =~ tr/ATGCatgc/TACGtacg/;
	return($rev);
}

sub reverseCigar {

	#invert cigar string so that corresponds to opposite strand
	my ($oriCig) = @_;

	#get and reverse letters and numbers from cigar
	
	my (@letters, @numbers);
	getCigarParts($oriCig, \@letters, \@numbers);
	@letters = reverse(@letters);
	@numbers = reverse(@numbers);
	
	#check that letters and numbers have same number of elements
	unless ((scalar @letters) == (scalar @numbers)) { print "could not reverse cigar ${oriCig}!\n"; return;}
	
	#reconstruct reversed cigar
	my $newCig;
	foreach my $number (@numbers) {
		$newCig .= ($number.shift(@letters));	
	}
	
	return $newCig;

}

sub zeroPad {
	#if we have a 1D array of start and stop locations of alignments, zero-pad these in order to correctly sort by start position
	#assume that 
	
	my @aligns = @_;
	
	#get number of digits in the longest start or stop coordinate
	#this is the number of total digits in all coordinates after zeropadding
	my $longest = 0;
	foreach my $element (@aligns) {
		
		#get start and end
		my ($start, $end, @other) = (split('xxx', $element));
		
		#check if number of digits in $start and $end is more than $longest
		if (length($start) > $longest) { $longest = length($start); }
		elsif (length($end) > $longest) { $longest = length($end); }
		
	}
	
	#add padding
	my @padded;
	foreach my $element (@aligns) {
		
		#get start and end
		my ($start, $end, @other) = (split('xxx', $element));
		
		$start = "0" x ($longest - length($start)) . $start;
		
		$end = "0" x ($longest - length($end)) . $end;
		
		push(@padded, join('xxx', $start, $end, @other));
	}
	
	return @padded;


}


1;
