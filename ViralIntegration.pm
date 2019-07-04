#!/usr/bin/perl

package ViralIntegration;
use Exporter;
@ISA = ('Exporter');
@EXPORT = qw(analyzeRead processCIGAR extractSeqCoords extractOutput extractCoords extractSequence getMatchedRegion getMatchedRegions getSecSup isAmbigLoc isRearrange printOutput printBed printMerged reverseComp reverseCigar);


##### subroutines #####

sub analyzeRead {
### Return the positiions flanking the integration site
### relative to the reference and the read
### returns:
###		intStart - integration start (relative to ref sequence)
###		intStop  - integration stop (relative to ref sequence)
###		relStart - integration start (relative to read sequence)
###		relStop  - integration stop (relative to read sequence)
	my ($start, $cigar, $end) = @_;
	
	my ($clipped,$mapped);
	if ($end eq "-") { ($clipped,$mapped) = ($cigar =~ /^(\d+)[SH](\d+)M/); }
	else 		 { ($mapped,$clipped) = ($cigar =~ /(\d+)M(\d+)[SH]$/); }

	my ($intStart, $intStop, $relStart, $relStop);
	if ($end eq "+") { 
		$intStart = $start + $mapped; 
		$intStop  = $start + $mapped + 1;
		$relStart = $mapped;
		$relStop  = $mapped + 1;
	}
	else {
		$intStop  = $start;
		$intStart = $start - 1;
		$relStop  = $clipped;
		$relStart = $clipped - 1;
	}
	return ($intStart, $intStop, $relStart, $relStop, $end)
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
	# This will give the toal amount of the read that was aligned (matched) and the part of the read that wasn't (clipped)
	# The assumption is that everything else is either matched or inserted (and therefore should still be counted)

	my $matched = length($seq) - $clipped;		
	
	# Construct simplified CIGAR string (only contains matched and clipped)
	if ($order == 1) { return(join('',($clipped,"S",$matched,"M"))); }
	else		 { return(join('',($matched,"M",$clipped,"S"))); }
}

sub extractSeqCoords {
### Extract coordinates of the human/viral sequences relative to read

	my ($ori, $dir, $alig, $overlap, $length) = @_;

	if ($dir eq "f") { # Read was aligned in fwd direction
		if ($ori eq "+") { return(1,($alig - $overlap)); }
		else		 { return(($length - $alig + $overlap + 1),$length); }
	} 
	else { # Read is oriented in the rev direction
		if ($ori eq "+") { return(($length - $alig + $overlap + 1),$length); }
		else	         { return(1,($alig - $overlap)); }
	}
}

sub extractOutput {
### Construct the output line
### Output line is a tab seperated string with the following information
### humanChr humanStart humanStop viralChr viralStar viralStop AmbiguousBases OverlapType Orientation(human) humanSeq viralSeq AmbiguousSeq HumanSecAlign ViralSecAlign
	my ($viralData) = shift @_;
	my ($humanData) = shift @_;
	my @intData     = @_;
	
	my @viral = split("\t",$viralData);
	my @human = split("\t",$humanData);

	my $vSec = $viral[-2];
	my $hSec = $human[-2];
	
	my $vSup = $viral[-1];
	my $hSup = $human[-1];

	my $overlap = $intData[6];
	
	my $overlaptype = $intData[8];
	
	my $vRearrange = $intData[9];
	my $hRearrange = $intData[10];
	
	my $hAmbig = $intData[11];
	my $vAmbig = $intData[12];

	### Extract junction coordinates relative to the target sequence
	my ($viralStart, $viralStop) = extractCoords($viral[1], $viral[2], $viral[5], $overlap);
	my ($humanStart, $humanStop) = extractCoords($human[1], $human[2], $human[5], $overlap);

	### Extract viral and human sequence componenets of the read
	### Account for overlap, this won't be included in either the human or viral segments
	my $viralSeq = substr($viral[6], $intData[4]-1, ($intData[5]-$intData[4]+1));
	my $humanSeq = substr($human[6], $intData[2]-1, ($intData[3]-$intData[2]+1));
	#my $viralSeq = extractSequence($intData[4], $intData[5]);
	#my $humanSeq = extractSequence($intData[2], $intData[3]);

	my $overlapSeq;
	if ($overlap > 0) { # only extract overlap sequence if there is one
		if ($intData[7] eq "hv") { $overlapSeq = substr($viral[6], $intData[3], $overlap); } # order along read is Human -> Viral
		else 			 { $overlapSeq = substr($viral[6], $intData[5], $overlap); } # order along read is Viral -> Human
	}
	else { $overlapSeq = ''; }# leave overlapSeq undefined


	#if ($viral[5] eq "+") 	{ $overlapSeq = substr($viral[6], ($viral[3]-$overlap), $overlap); }
	#else		  	{ $overlapSeq = substr($human[6], ($human[3]-$overlap), $overlap); }

	my $outline = join("\t", ($human[0], $humanStart, $humanStop, $viral[0], 
				  $viralStart, $viralStop, $overlap, $overlaptype, $human[5], 
				  $humanSeq, $viralSeq, $overlapSeq, $hSec, $vSec, $hRearrange, $vRearrange, $hAmbig, $vAmbig));

	return($outline);
}

sub extractCoords {
### Extract genomic coordinates for matched region based on start and stop relative to read
	my ($start, $stop, $ori, $overlap) = @_;

	if ($ori eq "-") { return(($start), ($stop + $overlap)); }
	else 		 { return(($start - $overlap), ($stop)); }	
}

sub extractSequence {
### Extract relevant sequence
	my ($seq, $start, $stop, $ori, $overlap) = @_;

	if ($ori eq "-") { return($seq = substr($seq, ($stop + $overlap))); }
	else			 { return($seq = substr($seq, 0, ($start - $overlap + 1))); }

}

sub getGenomicCoords {
	#get genomic coordinates of matched region (specified by cigar string, $cigar, and mapping orientation, $ori)
	#use left-most mapped position from sam file
	
	my ($pos, $cigar, $ori) = @_;
	
	#get matched regions from cigar
	my @aligns = getMatchedRegions($cigar, $ori);
	
	
	#need to know what the coordinate of the start of the read is, given the left-most mapping position
	
	
	my @coords;
	foreach my $align (@aligns) {
		my ($start, $stop) = split("xxx", $align);
	
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
	#returns a 1D array of matched regions in the form startxxxend	

	my ($cigar, $dir) = @_;
	if (($dir eq '-') or ($dir eq 'r')) { $cigar = reverseCigar($cigar); }
	
	#check for no alignments
	unless ($cigar =~ /M/) { return; }
	
	#if whole cigar is matched
	if ($cigar =~ /^(\d+)M$/) { my $match = ($cigar =~ /^(\d+)M$/); return "1xxx${match}"; }
	
	#get and reverse letters and numbers from cigar
	my @letters = split(/\d+/, $cigar);
	my @numbers = split(/[A-Z]/, $cigar);
	#since numbers precede letters, always get one extra empty element in letter array (which is at end)
	shift @letters;
	
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

sub isAmbigLoc {
### check for ambiguity in location of integration site in human or viral genome
## look to see if there are any secondary alignments that are equivalent to primary alignment

	#$dir is direction of primary alignment ('f' or 'r')
	#$cig is cigar of primary alignment
	#$sec is secondary alignments in the form /(chr,pos,CIGAR,NM;)*/
	#$aligL is number of bases mapped in primary alignment
	my ($dir, $cig, $sec) = @_;
	
	#get if match is at start or end of read
	
	
	#if alignment is in reverse orientation, convert to forward by reversing the cigar
	if ($dir eq 'r') { $cig = reverseCigar($cig); }
	
	#get if matched region is at start or end
	my ($primEnd, $aligL);
	if    (($cig =~ /^(\d+)[M]/)) { $aligL = $1; $primEnd = 'start'; }
	elsif (($cig =~ /(\d+)[M]$/)) { $aligL = $1; $primEnd = 'end';   }
	else 	{print "Can't figure out which end of the read is clipped in the primary alignment"; } #this shouldn't happen
	
	my ($secCigar, $secSense, $secEnd, $secAligL, $secAlignment);
	my $isAmbiguous = "";
	my @secs = split(";", $sec); # split secondary alignments to check one by one
	foreach $secAlignment (@secs) {
	
		if ($secAlignment eq "NA") { next; } #check for no secondary alignments
		my ($secSense, $secCigar) = (split(",",$secAlignment))[2,3]; #get sense and cigar from 
		
		if ($secCigar =~ /(^\d+[SH].*\d+[SH]$)/) { next; } #skip alignments where both ends are clipped
		if ($secCigar =~ /(^\d+[M].*\d+[M]$)/) { next; } #skip alignments where both ends are matched
		
		unless ($secCigar =~ /^\d+[M]|\d+[M]$/) { next; } #make sure one end is matched
	
		#check for if mapped is at start AND end (below assumes just one end mapped)
		#get if mapped part of supplementary alignment is at beginning or end of read
		if ($secSense eq '-') { $secCigar = reverseCigar($secCigar); } #reverse read if necessary
		
		if	  (($secCigar =~ /^(\d+)M/)) { $secEnd = 'start';  ($secAligL) = ($secCigar =~ /^(\d+)M/);}
		elsif (($secCigar =~ /(\d+)M$/)) { $secEnd = 'end'	 ;  ($secAligL) = ($secCigar =~ /[A-Z](\d+)M$/);}
		else { next; } #if mapped part not at start or end
		
		#if end of the read is the same and number of bases matched is the same, there is ambiguity
		if (($primEnd eq $secEnd) and ($secAligL == $aligL)) { $isAmbiguous = "yes"; } 
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
	
	my ($pCig, $pDir, $sup, $seq, $thresh) = @_;
	# $pCig is primary alignment of the first reference
	# $pDir is the direction of the primary alignment 
	# $sup are the secondary and supplementary alignments
	# $thresh is the fraction of the read that must be accounted for by alignments to be considered a rearrangement
	
	#if no sups or secs, can't be rearrangement
	if ($sup =~ /NA;NA/) { return ("no", 0, 0); }
	
	#if fully matched, can't be rearrangement
	if ($pCig =~ /^\d+M$/) { return ("no", 0, 0); }
	
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
	push(@aligns, getMatchedRegions($pCig, $pDir));
	
	#then do rest of the alignments
	my @supAligns = split(";", $sup);
	foreach my $supAlign (@supAligns) {
	
		if ($supAlign eq "NA") { next; } #check for no secondary alignments
		my ($supSense, $supCig) = (split(",",$supAlign))[2,3]; #get sense and cigar from alignment
		
		push(@aligns, getMatchedRegions($supCig, $supSense))
	}
	
	#sort array
	my @sorted = sort(zeroPad(@aligns));
	
	
	#remove any nested alignments
	#these are alignments that are completely encompassed by another alignment 
	#ie 001xxx100 and 010xxx020 are nested and 010xxx020 should be removed
	#also 001xxx009 and 001xxx010 are nested and 001xxx009 should be removed
	#keep any equivalent alignments - these are pairs with identical start and stop postions
	
	#strategy: start at end of array and check if the current alignment and the previous alignment are nested
	#if so, remove the nested one (which will always be the current alignment, since the list is sorted)
	my $lastKept = $#sorted;
	my @equivalent;
	for my $i (reverse((1..$#sorted))) {
		#get current and previous start and end
		my ($cS, $cE) = split('xxx', $sorted[$i]);
		my ($pS, $pE) = split('xxx', $sorted[$i-1]);
		#if two (or more) alignments are equivalent, we don't yet know if they should both/all be removed or not
		#just keep track of them
		if (($cS == $pS) and ($cE == $pE)) {
			push(@equivalent, $i);
			next;
		}
		#if current alignment and previous alignment are nested
		if ((($cS > $pS) and ($cE < $pE)) or (($cS > $pS) and ($cE == $pE))) { 
			#if no equivalent alignments
			if ($#equivalent == 0) {
				#just remove the one alignment (the current alignment)
				splice(@sorted, $i, 1); 
				}
			#if there are equivalent alignments
			else {
				#remove them all
				splice(@sorted, $i, ($#equivalent+1));
			}
		}
		#if current alignment and previous alignment are nested, and the previous alignment is the nested alignment
		#eg current = 001xxx120, previous = 001xxx100
		if (($cS == $pS) and ($cE > $pE)) {
			#first move current alignment to the position of the previous alignment, then delete
			@sorted[$i-1] = @sorted[$i];
			splice(@sorted, $i, 1);
		}
		#reset equivalent 
		@equivalent = ();
		
	}
	
	#check if alignments can account for whole read
	#check all alignments for gaps
	#to store the number of total bases gapped and number of gaps
	my ($gapBP, $gaps);
	$gapBP = (split('xxx', $sorted[0]))[0] - 1; #initialise to gap between start of read and start of first alignment
	if ($gapBP == 0) { $gaps = 0; }
	else { $gaps = 1; }
	
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
	
	my $isRearrange; #to store result, "yes" or "no"
	
	if ($gapBP < ($readlen - ($thresh * $readlen))) { $isRearrange = "yes"; }
	else { $isRearrange = "no"; }
	
	
	#get start and stop 
	
	return ($isRearrange, $gapBP, $gaps);
	
	
}

sub printOutput {
### Print output file
	my ($outFile, @outLines) = @_;	
	open (OUTFILE, ">$outFile") || die "Could not open output file: $outFile\n";
	print OUTFILE "Chr\tIntStart\tIntStop\tVirusRef\tVirusStart\tVirusStop\tNoAmbiguousBases\tOverlapType\tOrientation\tHostSeq\tViralSeq\tAmbiguousSeq\tHostSecondaryAlignments\tViralSecondaryAlignments\tPossibleHostTranslocation\tPossibleVectorRearrangement\tHostPossibleAmbiguous\tViralPossibleAmbiguous\tReadID\tmerged\n";
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
	
	my $temp = join("", ("temp.",$uni, ".bed"));
	open (TEMP, ">$temp");
	print TEMP $sorted;
	close TEMP;

	my $mergeOut = `bedtools merge -c 5,4 -o count,distinct -i $temp`;
	open (MERGED, ">$merged") || die "Could not open merged file: $merged\n";
	print MERGED $mergeOut;
	close MERGED;
	
	unlink($temp);
}

sub reverseComp {
### return reverse compliment of sequence
	my ($seq) = @_;
	my $rev = reverse $seq;
	$rev =~ tr/ATGCatgc/TACGtacg/;
	return($rev);
}

sub reverseCigar {

	my ($oriCig) = @_;

	#invert cigar string so that corresponds to opposite strand

	my (@letters, @numbers, $newCig, $letter, $number);

	#get and reverse letters and numbers from cigar
	@letters = reverse(split(/\d+/, $oriCig));
	@numbers = reverse(split(/[A-Z]/, $oriCig));
	#since numbers precede letters, always get one extra empty element in letter array (which is at end)
	pop @letters;
	
	#check that letters and numbers have same number of elements
	unless ((scalar @letters) == (scalar @numbers)) { print "could not reverse cigar ${oriCig}!\n"; return;}
	
	#reconstruct reversed cigar
	foreach $number (@numbers) {
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
		my ($start, $end) = (split('xxx', $element));
		
		#check if number of digits in $start and $end is more than $longest
		if (length($start) > $longest) { $longest = length($start); }
		elsif (length($end) > $longest) { $longest = length($end); }
		
	}
	
	#add padding
	my @padded;
	foreach my $element (@aligns) {
		
		#get start and end
		my ($start, $end) = (split('xxx', $element));
		
		$start = "0" x ($longest - length($start)) . $start;
		
		$end = "0" x ($longest - length($end)) . $end;
		
		push(@padded, join('xxx', $start, $end));
	}
	
	return @padded;


}


1;
