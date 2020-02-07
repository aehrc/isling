#!/usr/bin/perl
# It's not a bug, it's a feature

use strict;
use warnings;

use lib '.';

use ViralIntegration;
use Getopt::Long;



my $cutoff = 20; # each alignment must contain this number of aligned bases and this number of soft-clipped bases 
my $thresh = 0.95; #default amount of read that must be covered by alignments for rearrangement
my $tol = 1; #when processing CIGARS, combine any SHp IDPN elements between M regions with this number of bases or less
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
###	value:	{TargetID}xxx{ViralIntegrationData}xxx{Sequence}xxx{SeqOrientation}xxx{CIGAR}xxx{XA}xxx{SA}XX{NM}
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
	
	#get cigar
	my ($cig, $editDist2) = processCIGAR2($parts[5], $tol); # Note that could be a cigar or * if unmapped
	
	if ($cig =~ /(^\d+[SH].*\d+[SH]$)/) { next; } # skip alignments where both ends are clipped
	if ($cig =~ /[INDP]/) { next; } #don't consider any alignments with INDP elements
	unless ($cig =~ /^\d+[SH]|\d+[SH]$/) { next; } # at least one end of the read must be clipped in order to be considered
	unless ($cig) { next; } # keep checking to make sure double clipped reads don't sneak through
	
	#get supplementary (SA) and secondary (XA) alignments in order to check for possible vector rearrangements
	my ($vSec, $vSup) = getSecSup($vl);
	my $editDist = getEditDist($vl) + $editDist2;
	
	# if read is mapped in reverse orientation, reverse-complement the sequence
	my $seq;
	my $dir;
	if ($parts[1] & 0x10) {($seq) = reverseComp($parts[9]); $dir = 'r'; }
	else { $seq = $parts[9]; $dir = 'f'; }
	
	#get readID
	my $readID = $parts[0];
	
	#get 1-based mapping position:
	#from SAM format: 1-based leftmost mapping POSition of the first CIGAR operation that “consumes” a reference
	#base (see table below).  The first base in a reference sequence has coordinate 1. POS is set as 0 for
	#an unmapped read without coordinate.  If POS is 0, no assumptions can be made about RNAME and CIGAR
	my $pos = $parts[3];
	
	$viralIntegrations{join("xxx", ($readID, $seq))} = join("\t",($parts[2], $pos, $dir, $cig, $vSec, $vSup, $editDist));
	
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
	unless ($parts[5]) { next; }
	
	unless ($parts[5] =~ /^\d+[SH]|\d+[SH]$/) { next; }

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

	if (exists $viralIntegrations{join("xxx",($parts[0],$seq))}) { # only consider reads that were tagged from the viral alignment, no need to consider excess reads
		
		#get cigar
		my ($cig, $editDist2) = processCIGAR2($parts[5], $tol); # Note that could be a cigar or * if unmapped
		
		if ($cig =~ /(^\d+[SH].*\d+[SH]$)/) { next; } #skip alignments where both ends are clipped
		if ($cig =~ /[INDP]/) { next; } #don't consider any alignments with INDP elements
		
		unless ($cig) { next; }
		
		#get secondary and supplementary alignments from the line
		my ($hSec, $hSup) = getSecSup($hl);
		
		#get the edit distance, and add to it the bases that were modified when processing the CIGAR
		my $editDist = getEditDist($hl) + $editDist2;
		
		#get 1-based mapping position
		my $pos = $parts[3];
		
		#get readID
		my $readID = $parts[0];

		$humanIntegrations{join("xxx",($readID,$seq))} = join("\t",($parts[2], $pos, $dir, $cig, $hSec, $hSup, $editDist)); 
	}
}
close HUMAN;

### Look for evidence of integration sites by identifying chimeric reads
### Need to compare the junction sites in the viral and human genomes to see if there is any overlap (ambiguous bases)
### Adjust integration coordinates appropriately
my @outLines;
if ($verbose) { print "Detecting chimeric reads...\n"; }
foreach my $key (keys %viralIntegrations) {
	if (exists $humanIntegrations{$key}) { # only consider reads that are flagged in both human and viral alignments
		#for each putative integration, check that soft-clipped regions are complementary and calculate things for output
		my $outLine = collectIntersect($viralIntegrations{$key}, $humanIntegrations{$key}, $key, $thresh);
		#if (@intData) { $outLine = extractOutput($viralIntegrations{$key}, $humanIntegrations{$key}, @intData); }
		if ($outLine) { push(@outLines, join("\t", ($outLine, (split("xxx",$key))[0],(split("xxx",$key))[1]))) };
	}
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

sub printHelp {
	print "Pipeline for detection of viral integration sites within a genome\n\n";
	print "Usage:\n";
	print "\tperl softClip.pl --viral <sam> --human <sam> --cutoff <n> --thresh <n> --output <out> --bed <bed> --help\n\n";
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

sub collectIntersect {
### Check if there is overlap between the integration sites
### in human and viral integration sites
### returns integration start/stop relative to read sequence

### BWA Alignments are 1-based
	my ($viralData, $humanData, $key, $thresh) = @_;

	my ($vRef, $vPos, $vDir, $vCig, $vSec, $vSup, $vNM)  = (split("\t",$viralData));
	my ($hRef, $hPos, $hDir, $hCig, $hSec, $hSup, $hNM)  = (split("\t",$humanData));
	
	my ($readID, $seq) = split("xxx", $key);
	
	## orientation of integration site: + is after aligned sequence, - is before aligned sequence
	### ie if CIGAR is mapped/clipped orientation is +
	### and if CIGAR is clipped/mapped orientation is -
	### but consider direction too:
	### CIGAR strings are always reported relative to the strand
	### 100M50S on a fwd read = 50S100M on a rev read
	my $vOri = getOri($vCig, $vDir);
	my $hOri = getOri($hCig, $hDir);
	
	unless 	($vOri and $hOri) { return; } # Catch a weird case where an orientation isn't found. Appears to happen when both ends are clipped
	
	#if both reads have same part of read clipped/mapped, this can't be an integration
	if	($vOri eq $hOri) { return; } 

	### First find if there is any overlap between the human and viral junctions
	### Do so by comparing the CIGAR strings
	### Can calculate by subtracting aligned of one from clipped of other	
	my ($hClip) = ($hCig =~ /(\d+)S/);
	my ($hAlig) = ($hCig =~ /(\d+)M/);
	my ($vClip) = ($vCig =~ /(\d+)S/);
	my ($vAlig) = ($vCig =~ /(\d+)M/);

	#check that human and viral soft-clipped and aligned regions meet cutoff
	unless (($hClip > $cutoff) and ($hAlig > $cutoff) and ($vClip > $cutoff) and ($vAlig > $cutoff)) { return; }
	
	### Overlap should be the same regardless of how it's calculated so double check
	my $overlap1 = abs($hAlig - $vClip);
	my $overlap2 = abs($vAlig - $hClip);

	#here check that absolute values of overlap are the same
	unless (abs($overlap1) == abs($overlap2)) { 
		print "Impossible overlap found\n";
		$DB::single = 1;
		return;
	}
	my $overlap = abs($overlap1);
	
	#ambigous bases may result from either an overlap of aligned regions, or a gap
	#if the total number of aligned bases (from human and viral alignments) are greater than the read length, then it's an overlap
	#otherwise it's a gap
	my ($overlaptype, $readlen);
	$readlen = length($seq);
	if 		(($hAlig + $vAlig) > ($readlen))  {	$overlaptype = "overlap";	} #overlap
	elsif 	(($hAlig + $vAlig) == ($readlen)) {	$overlaptype = "none";		} #no gap or overlap
	else 									  {	$overlaptype = "gap";		} #gap

	#also enforce a cutoff on the number of unambiguously mapped bases
	#that is, if there is an overlap, the number of mapped bases excluding the overlapped region must still
	#be more than the cutoff
	if ($overlaptype eq "gap") {
		unless ((($vAlig - $overlap) > $cutoff) and (($hAlig - $overlap) > $cutoff)) { return; }		
	}
	
	#caculate total edit distance - sum of host and virus edit distance, and gap if there is one
	my $totalNM = $hNM + $vNM;
	if ($overlaptype eq 'gap') { $totalNM += $overlap; }
	
	#check to see is whole read can be accounted for by alignments to the vector
	my ($isVecRearrange);
	#decide if read is more likely a rearrangement or integration based on edit distance
	#if edit distances are the same, err on the side of caution and assign it as a rearrangement
	$isVecRearrange = isRearrangeOrInt($vCig, $vDir, $vRef, $vPos, $vSec, $vSup, $seq, $thresh, $vNM, $totalNM);

	my ($isHumRearrange);
	$isHumRearrange = isRearrangeOrInt($hCig, $hDir, $hRef, $hPos, $hSec, $hSup, $seq, $thresh, $hNM, $totalNM);
	
	#check to see if location of human alignment is ambiguous: multiple equivalent alignments accounting for human part of read
	my $isHumAmbig;
	if ($hSec eq "NA") { $isHumAmbig = "no";}
	else { $isHumAmbig = isAmbigLoc($hDir, $hCig, $hSec, 'soft', $seq, $tol);}
	
	#check to see if location of viral alignment is ambiguous: multiple equivalent alignments accounting for viral part of read
	my $isVirAmbig;
	if ($vSec eq "NA") { $isVirAmbig = "no";}
	else { $isVirAmbig = isAmbigLoc($vDir, $vCig, $vSec, 'soft', $seq, $tol);}

	### Calculate the start and stop positions of the viral and human sequences relative to the read
	my ($hRStart,$hRStop) = extractSeqCoords($hOri, $hAlig, $overlap, $readlen, $overlaptype);
	my ($vRStart,$vRStop) = extractSeqCoords($vOri, $vAlig, $overlap, $readlen, $overlaptype);
	
	### Collect integration start/stop positions relative to read
	### These are the bases that flank the bond that is broken by insertion of the virus
	### Takes any overlap into account	
	my ($intRStart, $intRStop, $order);

	if    ($hOri eq "+") { ($intRStart,$intRStop,$order) = ($hRStop+1,$vRStart,"hv"); } # Orientation is Human -> Virus
	elsif ($hOri eq "-") { ($intRStart,$intRStop,$order) = ($vRStop+1,$hRStart,"vh"); } # Orientation is Virus -> Human
	else 			   { 
		print "Something weird has happened";
		return;
	}
	
	### Extract viral and human sequence componenets of the read
	### Account for overlap, this won't be included in either the human or viral segments
	
	my $viralSeq = substr($seq, $vRStart-1, $vRStop-$vRStart+1);
	my $humanSeq = substr($seq, $hRStart-1, $hRStop-$hRStart+1);
	my $overlapSeq = substr($seq, $intRStart, $overlap);	
	
	#check that we have extracted the right number of bases (sum of three regions should be whole read)
	if ((length($viralSeq) + length($humanSeq) + length($overlapSeq)) != $readlen) { print "wrong number of bases extracted\n"; }

	### Extract junction coordinates relative to the target sequence
	
	#first get start and stop genomic coordinates of integration
	my ($hgStart, $hgStop) = extractCoords($hAlig, $overlap, $overlaptype, $hPos, $hOri);
	my ($vgStart, $vgStop) = extractCoords($vAlig, $overlap, $overlaptype, $vPos, $vOri);
	

	#generate output
	my $outline = join("\t", ($hRef, $hgStart, $hgStop, $vRef, 
				  $vgStart, $vgStop, $overlap, $overlaptype, $order, 
				  $humanSeq, $viralSeq, $overlapSeq, $hNM, $vNM, $totalNM, $isHumRearrange, 
				  $isVecRearrange, $isHumAmbig, $isVirAmbig, 'chimeric'));

}

sub extractOutput {
### Construct the output line
### Output line is a tab seperated string with the following information
### humanChr humanStart humanStop viralChr viralStar viralStop AmbiguousBases OverlapType Orientation(human) humanSeq viralSeq AmbiguousSeq HumanSecAlign ViralSecAlign
	my ($viralData) = shift @_;
	my ($humanData) = shift @_;
	my @intData     = @_;
	
	#viralData and $humanData have the following fields:
	#$refName, $pos, $dir, $cig, $vSec, $vSup, $editDist
	
	#intData has the fields
	#$intRStart, $intRStop, $hRStart, $hRStop, $vRStart, $vRStop, $overlap1, $order, $overlaptype, $isVecRearrange, $isHumRearrange, $isHumAmbig, $isVirAmbig
	
	my @viral = split("\t",$viralData);
	my @human = split("\t",$humanData);

	my $overlap = $intData[6];
	my $overlaptype = $intData[8];
	
	my $vRearrange = $intData[9];
	my $hRearrange = $intData[10];
	
	my $hAmbig = $intData[11];
	my $vAmbig = $intData[12];
	
	my $vNM = $viral[-1];
	my $hNM = $human[-1];
	
	my $totalNM = $hNM + $vNM;
	if ($overlaptype eq 'gap') { $totalNM += $overlap; }

	### Extract junction coordinates relative to the target sequence
	my ($viralStart, $viralStop) = extractCoords($viral[1], $viral[2], $viral[5], $overlap);
	my ($humanStart, $humanStop) = extractCoords($human[1], $human[2], $human[5], $overlap);

	### Extract viral and human sequence componenets of the read
	### Account for overlap, this won't be included in either the human or viral segments
	my $viralSeq = substr($viral[6], $intData[4]-1, ($intData[5]-$intData[4]+1));
	my $humanSeq = substr($human[6], $intData[2]-1, ($intData[3]-$intData[2]+1));

	my $overlapSeq;
	if ($overlap > 0) { # only extract overlap sequence if there is one
		if ($intData[7] eq "hv") { $overlapSeq = substr($viral[6], $intData[3], $overlap); } # order along read is Human -> Viral
		else 			 { $overlapSeq = substr($viral[6], $intData[5], $overlap); } # order along read is Viral -> Human
	}
	else { $overlapSeq = ''; }# leave overlapSeq undefined


	my $outline = join("\t", ($human[0], $humanStart, $humanStop, $viral[0], 
				  $viralStart, $viralStop, $overlap, $overlaptype, $intData[7], 
				  $humanSeq, $viralSeq, $overlapSeq, $hNM, $vNM, $totalNM, $hRearrange, $vRearrange, $hAmbig, $vAmbig, 'chimeric'));

	return($outline);
}


sub getOri {
	#get orientation of soft-clipped alignment from CIGAR
	
	## orientation = orientation of integration site: + is after aligned sequence, - is before aligned sequence
	### ie if CIGAR is mapped/clipped orientation is +
	### and if CIGAR is clipped/mapped orientation is -
	
	#note that CIGAR is reversed if read is mapped in reverse orientation, so 50M20S if forward is equivalent to 20S50M if reverse
	
	my ($cig, $dir) = @_;
	
	my $ori;
	
	#if read is forward
	if ($dir eq 'f') {
		if ($cig =~ /(^\d+)[SH]/)   	  { $ori = "-"; }  # integration site is before the viral sequence
		elsif ($cig =~ /(\d+)[SH]$/) 	  { $ori = "+"; } # integration site is after the viral sequence
	}
	
	#if read is reverse
	elsif ($dir eq 'r') {
		if ($cig =~ /(^\d+)[SH]/)   	  { $ori = "+"; }  # integration site is before the viral sequence
		elsif ($cig =~ /(\d+)[SH]$/) 	  { $ori = "-"; } # integration site is after the viral sequence

	}
	
	return $ori;
}
