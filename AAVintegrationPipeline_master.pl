#!/usr/bin/perl
# It's not a bug, it's a feature

use strict;
use warnings;

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
### These hashe's contain the information about the integration sites within the resepctive genomes
### Data is saved as a tab seperated string
### Read name => Chr intStart intStop relStart relStop orientation readSeqeunce


### Collect junction reads from viral genome
### Junctions are defined by the CIGAR string: SM or MS
### Save the resulting information in a hash:
###	key: 	{ReadID}xxx{ReadSequence}
###	value:	{TargetID}xxx{ViralIntegrationData}xxx{Sequence}xxx{SeqOrientation}xxx{CIGAR}xxx{XA}
###		SeqOrientation: f = aligned in fwd orientation
###				r = aligned in rev orientation
###		Sequence is always given in fwd orientation
### XA is the secondary alignments in the XA field (from BWA); not present use "noXA"

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

	#Add get secondary alignment information from XA field
	my $vSec;
	if ($vl =~ /XA\:.+\:.+/) 	{ 
		($vSec) = ($vl =~ /XA\:.+\:(.+)\s/); 
	}
	else 						{ $vSec = "NA"; }
	
	if   (@viralInt and ($parts[1] & 0x10)) { $viralIntegrations{join("xxx", ($parts[0],(reverseComp($parts[9]))[0]))} = join("\t",($parts[2], @viralInt, (reverseComp($parts[9]))[0], 'r', $cig, $vSec)); }
	elsif (@viralInt) 			{ $viralIntegrations{join("xxx", ($parts[0],$parts[9]))}	           = join("\t",($parts[2], @viralInt, $parts[9],                   'f', $cig, $vSec)); } 
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
	
		#get secondary alignments from XA field
		my $hSec;
		if ($hl =~ /XA\:.+\:.+/) 	{ ($hSec) = ($hl =~ /XA\:.+\:(.+)\s/); }
		else 						{ $hSec = "NA"; }
		
		my @humanInt;
		if    ($cig =~ /^\d+[SH].+\d+[SH]$/) { @humanInt = undef; }
		elsif ($cig =~ /(^\d+)[SH]/)         { if ($1 > $cutoff) { @humanInt = analyzeRead($parts[3], $cig, "-"); } }
		elsif ($cig =~ /(\d+)[SH]$/)         { if ($1 > $cutoff) { @humanInt = analyzeRead($parts[3], $cig, "+"); } }

		if (@humanInt) { $humanIntegrations{join("xxx",($parts[0],$seq))} = join("\t",($parts[2], @humanInt, $seq, $ori, $cig, $hSec)); }
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
	elsif ($oriCig =~ /^(.+)(\d+)[SH]$/) { ($other, $clipped, $order) = ($1,$2,2); }
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

	my ($vStart, $vStop, $vOri, $seq, $vDir, $vCig) = (split("\t",$viralData))[3,4,5,6,7,8];
	my ($hStart, $hStop, $hOri, $hDir, $hCig)       = (split("\t",$humanData))[3,4,5,7,8];

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
	my $overlaptype;
	if 		(($hAlig + $vAlig) > ($hClip + $hAlig))	 {	$overlaptype = "overlap";	} #overlap
	elsif 	(($hAlig + $vAlig) == ($hClip + $hAlig)) {	$overlaptype = "none";		} #no gap or overlap
	else 											 {	$overlaptype = "gap";		} #gap
	


	### Calculate the start and stop positions of the viral and human sequences relative to the read
	my ($hRStart,$hRStop) = extractSeqCoords($hOri, $hDir, $hAlig, abs($overlap1), length($seq));
	my ($vRStart,$vRStop) = extractSeqCoords($vOri, $vDir, $vAlig, abs($overlap1), length($seq));

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

	return($intRStart, $intRStop, $hRStart, $hRStop, $vRStart, $vRStop, $overlap1, $order, $overlaptype);

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

	my $vSec = pop(@viral);
	my $hSec = pop(@human);

	my $overlap = $intData[6];
	
	my $overlaptype = $intData[8];

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
				  $humanSeq, $viralSeq, $overlapSeq, $hSec, $vSec));

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

sub printOutput {
### Print output file
	my ($outFile, @outLines) = @_;	
	open (OUTFILE, ">$outFile") || die "Could not open output file: $outFile\n";
	print OUTFILE "Chr\tIntStart\tIntStop\tVirusRef\tVirusStart\tVirusStop\tNo. Ambiguous Bases\tOverlap Type\tOrientation\tHuman Seq\tViral Seq\tAmbiguous Seq\tHuman Secondary Alignments\tViral Secondary Alignments\tRead ID\t(merged)\n";
	foreach my $line (@outLines) { print OUTFILE "$line\n"; }
	close OUTFILE;
}

sub printBed {
### Print Bed output
### Bed starts are 0-based and there ends are 1-based
### 1	0	1 <- referes to the first base on chromosome 1

	my ($bedFile, @outLines) = @_;
	open (BEDFILE, ">$bedFile") || die "Could not open bedfile: $bedFile\n";
	foreach my $line (@outLines) {
		my @parts = split("\t",$line);
		print BEDFILE "$parts[0]\t$parts[1]\t$parts[2]\t$parts[14]\t.\t$parts[8]\n";
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
