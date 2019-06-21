#!/usr/bin/perl

package ViralIntegration;
use Exporter;
@ISA = ('Exporter');
@EXPORT = qw(analyseRead processCIGAR collectIntersect extractSeqCoords extractOutput extractCoords extractSequence getSecSup isAmbigLoc isRearrange printOutput printBed printMerged reverseComp reverseCigar);


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

sub collectIntersect {
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
	#if 	($vOri eq $hOri)  { 
	#	return; 
	#} # skip reads where same side is clipped in virus and human (not a true junction)

	#if ($vDir eq 'r' and $hDir eq 'f' and $vOri eq "+") { 
	#	$vOri   = "-";
	#	$vStart = length($seq) - $vStart - 1;
	#	$vStop  = length($seq) - $vStop  - 1;
	#}
	#elsif ($vDir eq 'r' and $hDir eq 'f' and $vOri eq "-") { 
	#	$vOri   = "+";
	#	$vStart = length($seq) - $vStart - 1;
	#	$vStop  = length($seq) - $vStop  - 1;
	#}
	#elsif ($vDir eq 'f' and $hDir eq 'r' and $hOri eq "+") { 
	#	$hOri   = "-"; 
	#	$hStart = length($seq) - $hStart - 1;
	#	$hStop  = length($seq) - $hStop  - 1; 
	#}
	#elsif ($vDir eq 'f' and $hDir eq 'r' and $hOri eq "-") { 
	#	$hOri   = "+"; 
	#	$vStart = length($seq) - $vStart - 1;
	#	$vStop  = length($seq) - $vStop  - 1; 
	#}

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

	#my ($hRStart, $hRStop);
	#if ($hOri eq "+") { # human integration site is after the aligned sequence: Human -> Viral
	#	$hRStart = 1; # because the first half of the read is human
	#	$hRStop  = $hAlig - $overlap1; # it extends to the start of the overlap
	#}
	#else { # integration is before the aligned seqeunce: Viral -> Human
	#	$hRStop  = length($seq); # because the second half of the read is human
	#	$hRStart = $hRStop - $hAlig + $overlap1 + 1; #extends to the start of the overlap
	#}
		
	#repeat for the viral junction
	#my ($vRStart, $vRStop);
	#if ($vOri eq "+") {
	#	$vRStart = 1;
	#	$vRStop  = $vAlig - $overlap1;
	#}
	#else {
	#	$vRStop  = length($seq);
	#	$vRStart = $vRStop - $vAlig + $overlap1 + 1; 
	#}

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

	#if ($vOri eq "-" and $hOri eq "+") { # Orientation is Human -> Virus
	#	$intRStart = $hRStop;
	#	$intRStop  = $vRStart;
	#	$order     = "hv";
	#}
	#elsif ($hOri eq "-" and $vOri eq "+") { # Orientation is Virus -> Human
	#	$intRStart = $vRStop;
	#	$intRStop  = $vRStart;
	#	$order     = "vh";
	#}

	#else {
	#	return; # If you end up here something's gone weird. Need to double check to make sure nothing would end up here
	#}

	return($intRStart, $intRStop, $hRStart, $hRStop, $vRStart, $vRStop, $overlap1, $order, $overlaptype, $isVecRearrange, $isHumRearrange, $isHumAmbig, $isVirAmbig);

	#if 	  ($vOri eq "+" and $hOri eq "-" and $vStop - 1 > $hStart) { return($hStart, $vStop); } # overalp with order virus -> human
	#elsif ($vOri eq "-" and $hOri eq "+" and $hStop - 1 > $vStart)    { return($vStart, $hStop); } # overlap with order human -> virus
	#else 							   	   { return($hStart, $hStop); } # if no overlap, use human positions
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
### Extract coordinates
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

sub getSecSup {
#get secondary (XA) and supplementary (SA) alignments from line

#note that secondary alignemnts have the form /XA:Z:(chr,pos,CIGAR,NM;)*/
#supplementary alignments have the form /SA:Z:(chr,pos,strand,CIGAR,PHRED,NM;)*/

	my ($line) = @_;
	
	#get supplementary alignments from SA field
	my ($sup, $sec);
	if ($line =~ /SA:.:.*;/) 	{ 
		($sup) = ($line =~ /SA:.:(.*?);/); #will this get only one supp alignment?
	
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
	my ($dir, $cig, $sec, $aligL) = @_;
	
	#get if match is at start or end of read
	
	
	#if alignment is in reverse orientation, convert to forward by reversing the cigar
	if ($dir eq 'r') { $cig = reverseCigar($cig); }
	
	#get if matched region is at start or end
	my $primEnd;
	if    (($cig =~ /^(\d+)[M]/)) { $primEnd = 'start'; }
	elsif (($cig =~ /(\d+)[M]$/)) { $primEnd = 'end';   }
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

	my ($pCig, $pDir, $sup, $readlen) = @_;
	
	#after processing, cigars have one matched region per alignment
	#get location in read of matched region and its length
	#make array of start of matched regions and their lengths in format startxxxend
	#need to take into account direction of read - convert everything to forward orientation
	#need to zeropad start and length so that each have three digits, for sorting later
	#NOTE - THIS ASSUMES READ LENGTHS LESS THAN 1000!!!
	if ($readlen >= 1000) { print "warning: vector rearangement check doesn't work for read lengths >= 1000\n" ;}
	
	
	#note - currently considering all alignments together regardless of their reference: consider changing to check only that alignments to a sigle reference account for the whole read
	#for vector/virus we would want to do this because probably don't have different viruses recombining?
	#for human/host we would not want to do this could have translocations
	#could do this by making one @aligns per reference, and store all 
	#note - using 1-based indexing of read alignments
		
	my @aligns; #to store start and end of each alignment
	
	#first do primary alignment
	if ($pDir eq "r") { $pCig = reverseCigar($pCig); } #reverse cigar if necessary
	#primary alignment cigar has been processed so we can assume it only has one mapped region and one soft-clipped region
	if 	   ($pCig =~ /^(\d+)[M]/) { push(@aligns, "1xxx${1}"); }  #simplest case is that matched region is start
	elsif  ($pCig =~ /^(\d+)[SH](\d+)[M]$/) { push(@aligns, (${1}+1)."xxx".(${1}+${2})); }
	
	#now do rest of the alignments
	my @supAligns = split(";", $sup);
	foreach my $supAlign (@supAligns) {
	
		if ($supAlign eq "NA") { next; } #check for no secondary alignments
		my ($supSense, $supCig) = (split(",",$supAlign))[2,3]; #get sense and cigar from alignment
	
		if ($supSense eq '-') { $supCig = reverseCigar($supCig); }
		
		#get matched region(s) from alignment
		#get and reverse letters and numbers from cigar
		my @letters = split(/\d+/, $supCig);
		my @numbers = split(/[A-Z]/, $supCig);
		#since numbers precede letters, always get one extra empty element in letter array (which is at start)
		shift @letters;


		#for each match in @letters, get position in read where alignment starts and length of alignment
		my ($start, $matchedLen, $end);
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
	}
	
	#add padding and sort alignments
	#zeropadding is needed so that sorting works correctly
	my @padded;
	foreach my $element (@aligns) {
		
		#get start and end
		my ($start, $end) = (split('xxx', $element));
		
		#add padding to start if not three digits
		if ($start =~ /^(\d)$/) { $start = "00${1}"; }
		elsif ($start =~ /^(\d\d)$/) {$start = "0${1}"; }
		
		#add padding to end if not three digits
		if ($end =~ /^(\d)$/) { $end = "00${1}"; }
		elsif ($end =~ /^(\d\d)$/) {$end = "0${1}"; }
		
		push(@padded, join('xxx', $start, $end));
	}
	
	#sort array
	my @sorted = sort(@padded);
	
	#check if alignments can account for whole read
	my $isRearrange; #to store result, "yes" or "no"
	#check that first base is accounted for
	#since list is sorted, if first base is accounted for then first alignment will account for it
	if ($sorted[0] !~ "001") { return "no"; } 

	#check rest of the read
	my $currentPos = 0;
	my $lastUsed = -1;
	my $i = (scalar @sorted) - 1; #start at end of list and go backwards
	
	while (($currentPos < $readlen) and ($i != $lastUsed )) {
		if ((split('xxx', $sorted[$i]))[0] <= ($currentPos +1))  { 
			$currentPos = (split('xxx', $sorted[$i]))[1]; #update current position
			$lastUsed = $i; #update last used alignment
			$i = ((scalar @sorted) - 1); #start again at end of list	
		}
		else { $i -= 1; }
	}
	
	#check if we reached the end of the read
	$isRearrange = ($currentPos == $readlen) ? "yes" : "no";

	return $isRearrange;

}

sub printOutput {
### Print output file
	my ($outFile, @outLines) = @_;	
	open (OUTFILE, ">$outFile") || die "Could not open output file: $outFile\n";
	print OUTFILE "Chr\tIntStart\tIntStop\tVirusRef\tVirusStart\tVirusStop\tNoAmbiguousBases\tOverlapType\tOrientation\tHostSeq\tViralSeq\tAmbiguousSeq\tHostSecondaryAlignments\tViralSecondaryAlignments\tPossibleHostTranslocation\tPossibleVectorRearrangement\tHostPossibleAmbiguous\tViralPossibleAmbiguous\tReadID\tMergedRead\n";
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


1;
