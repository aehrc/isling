#!/usr/bin/perl
# It's not a bug, it's a feature

use strict;
use warnings;

use lib '.';

use ViralIntegration;
use Getopt::Long;



my $cutoff = 20; # each alignment must contain this number of aligned bases and this number of soft-clipped bases 
my $thresh = 0.95; #default amount of read that must be covered by alignments for rearrangement
my $tol = 5; #when processing CIGARS, combine any SHp IDPN elements between M regions with this number of bases or less
my $viral;
my $host;
my $output = "integrationSites.txt";
my $bed;
my $merged;
my $verbose;
my $help;

GetOptions('cutoff=i' => \$cutoff,
		   'thresh=f' => \$thresh,
		   'viral=s'  => \$viral,
		   'host=s'  => \$host,
		   'output=s' => \$output,
		   'bed=s'    => \$bed,
		   'merged=s' => \$merged,
		   'tol=i'    => \$tol,
		   'verbose'  => \$verbose,
		   'help',    => \$help);

if ($help) { printHelp(); }

unless ($viral and $host) { printHelp(); }

# check inputs
if ($cutoff < 0) { die "Cutoff must be greater than zero\n"; }
if ($tol < 0) { die "Tolerance must be greater than zero\n"; }

my %viralIntegrations;
my %hostIntegrations;
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
	my ($cig, $extraNM) = simplifyCIGAR($parts[5], $parts[9], $tol);
	
	unless ($cig) { next; } # keep checking to make sure double clipped reads don't sneak through
	if ($cig =~ /(^\d+[SH].*\d+[SH]$)/) { next; } # skip alignments where both ends are clipped
	if ($cig =~ /[INDP]/) { next; } #don't consider any alignments with INDP elements
	unless ($cig =~ /^\d+[SH]|\d+[SH]$/) { next; } # at least one end of the read must be clipped in order to be considered
	
	
	#get supplementary (SA) and secondary (XA) alignments in order to check for possible vector rearrangements
	my ($vSec, $vSup) = getSecSup($vl);
	my $editDist = getEditDist($vl) + $extraNM;
	
	# if read is mapped in reverse orientation, reverse-complement the sequence
	my $seq;
	my $dir;
	if ($parts[1] & 0x10) {($seq) = reverseComp($parts[9]); $dir = 'r'; }
	else { $seq = $parts[9]; $dir = 'f'; }
	
	#get readID
	my $readID = $parts[0];
	my $mapq = $parts[4];
	
	# append R1 or R2 to $ID if flag is set
	if ($parts[1] & 0x40) { $readID = $readID."/1"; }
	if ($parts[1] & 0x80) { $readID = $readID."/2"; }
	
	#get 1-based mapping position:
	#from SAM format: 1-based leftmost mapping POSition of the first CIGAR operation that “consumes” a reference
	#base (see table below).  The first base in a reference sequence has coordinate 1. POS is set as 0 for
	#an unmapped read without coordinate.  If POS is 0, no assumptions can be made about RNAME and CIGAR
	my $pos = $parts[3];
	
	$viralIntegrations{join("xxx", ($readID, $seq))} = join("\t",($parts[2], $pos, $dir, $cig, $vSec, $vSup, $editDist, $mapq));
	
}
close VIRAL;

### Collect junction reads from host genome
### This is the same process as for viral junctions
open (HOST, $host) || die "Could not open host alignment file: $host\n";
if ($verbose) { print "Processing host alignment...\n"; }
while (my $hl = <HOST>) {
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
	
	#get readID
	my $readID = $parts[0];
	
	# append R1 or R2 to $readID if flag is set
	if ($parts[1] & 0x40) { $readID = $readID."/1"; }
	if ($parts[1] & 0x80) { $readID = $readID."/2"; }

	if (exists $viralIntegrations{join("xxx",($readID,$seq))}) { # only consider reads that were tagged from the viral alignment, no need to consider excess reads
		
		#get cigar
		my ($cig, $extraNM) = simplifyCIGAR($parts[5], $parts[9], $tol);
		
		unless ($cig) { next; }
		if ($cig =~ /(^\d+[SH].*\d+[SH]$)/) { next; } #skip alignments where both ends are clipped
		if ($cig =~ /[INDP]/) { next; } #don't consider any alignments with INDP elements
		
		#get secondary and supplementary alignments from the line
		my ($hSec, $hSup) = getSecSup($hl);
		
		#get the edit distance, and add to it the bases that were modified when processing the CIGAR
		my $editDist = getEditDist($hl) + $extraNM;
		
		#get 1-based mapping position
		my $pos = $parts[3];
		
		my $mapq = $parts[4];

		$hostIntegrations{join("xxx",($readID,$seq))} = join("\t",($parts[2], $pos, $dir, $cig, $hSec, $hSup, $editDist, $mapq)); 
	}
}
close HOST;

### Look for evidence of integration sites by identifying chimeric reads
### Need to compare the junction sites in the viral and host genomes to see if there is any overlap (ambiguous bases)
### Adjust integration coordinates appropriately
my @outLines;
if ($verbose) { print "Detecting chimeric reads...\n"; }
foreach my $key (keys %viralIntegrations) {
	if (exists $hostIntegrations{$key}) { # only consider reads that are flagged in both host and viral alignments
		#for each putative integration, check that soft-clipped regions are complementary and calculate things for output
		my $outLine = collectIntersect($viralIntegrations{$key}, $hostIntegrations{$key}, $key, $thresh);
		#if (@intData) { $outLine = extractOutput($viralIntegrations{$key}, $hostIntegrations{$key}, @intData); }
		if ($outLine) { push(@outLines, join("\t", ($outLine, (split("xxx",$key))[0],(split("xxx",$key))[1]))) };
	}
}


if ($verbose) { print "Writing output...\n"; }

my $header = "Chr\tIntStart\tIntStop\tVirusRef\tVirusStart\tVirusStop\tNoAmbiguousBases\tOverlapType\tOrientation\tVirusOrientation\tHostSeq\tVirusSeq\tAmbiguousSeq\tHostEditDist\tViralEditDist\tTotalEditDist\tPossibleHostTranslocation\tPossibleVectorRearrangement\tHostPossibleAmbiguous\tViralPossibleAmbiguous\tType\tHostMapQ\tViralMapQ\tReadID\tReadSeq\n";		

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
	print "\tperl softClip.pl --viral <sam> --host <sam> --cutoff <n> --thresh <n> --tol <n> --min-mapq <n> --output <out> --bed <bed> --help\n\n";
	print "Arguments:\n";
	print "\t--viral:   Alignment of reads to viral genomes (sam)\n";
	print "\t--host:   Alignment of reads to host genome (sam)\n";
	print "\t--cutoff:  Minimum number of clipped reads to be considered (default = 20)\n";
	print "\t--thresh:	Amount of read that must be covered by one alignment to be considered rearrangement";
	print "\t--tol:	Soft-clipped regions at either end of the read shorter than this will be absorbed into the nearest mapped region";
	print "\t--min-mapq:	Minimum mapping quality to consider a read";	
	print "\t--output:  Output file for results (default = integrationSite.txt\n";
	print "\t--bed:     Print integrations sites to indicated bed file (default = NA)\n";
	print "\t--merged:  Merge bedfile into overlapping integration sites (default = NA)\n";
	print "\t--verbose: Print progress messages\n";
	print "\t--help:    Print this help message and quit\n";

	exit;
}

sub collectIntersect {
### Check if there is overlap between the integration sites
### in host and viral integration sites
### returns integration start/stop relative to read sequence

### BWA Alignments are 1-based
	my ($viralData, $hostData, $key, $thresh) = @_;

	my ($vRef, $vPos, $vDir, $vCig, $vSec, $vSup, $vNM, $vMapq)  = (split("\t",$viralData));
	my ($hRef, $hPos, $hDir, $hCig, $hSec, $hSup, $hNM, $hMapq)  = (split("\t",$hostData));
	
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
	
	# if both reads have same part of read clipped/mapped, this can't be an integration
	if	($vOri eq $hOri) { return; } 

	### First find if there is any overlap between the host and viral junctions
	### Do so by comparing the CIGAR strings
	### Can calculate by subtracting aligned of one from clipped of other	
	my ($hClip) = ($hCig =~ /(\d+)S/);
	my ($hAlig) = ($hCig =~ /(\d+)M/);
	my ($vClip) = ($vCig =~ /(\d+)S/);
	my ($vAlig) = ($vCig =~ /(\d+)M/);

	# check that host and viral aligned regions meet cutoff
	unless (($hAlig >= $cutoff) and ($vAlig >= $cutoff)) { return; }
	
	### Overlap should be the same regardless of how it's calculated so double check
	my $overlap1 = abs($hAlig - $vClip);
	my $overlap2 = abs($vAlig - $hClip);

	# here check that absolute values of overlap are the same
	unless (abs($overlap1) == abs($overlap2)) { 
		print "Impossible overlap found\n";
		return;
	}
	my $overlap = abs($overlap1);
	
	# ambigous bases may result from either an overlap of aligned regions, or a gap
	# if the total number of aligned bases (from host and viral alignments) are greater than the read length, then it's an overlap
	# otherwise it's a gap
	my ($overlaptype, $readlen);
	$readlen = length($seq);
	if 		(($hAlig + $vAlig) > ($readlen))  {	$overlaptype = "overlap";	} #overlap
	elsif 	(($hAlig + $vAlig) == ($readlen)) {	$overlaptype = "none";		} #no gap or overlap
	else 									  {	$overlaptype = "gap";		} #gap

	# also enforce a cutoff on the number of unambiguously mapped bases
	# that is, if there is an overlap, the number of mapped bases excluding the overlapped region must still
	# be more than the cutoff
	if ($overlaptype eq "overlap") {
		unless ((($vAlig - $overlap) >= $cutoff) and (($hAlig - $overlap) >= $cutoff)) { return; }		
	}
	
	# caculate total edit distance - sum of host and virus edit distance, and gap if there is one
	my $totalNM = $hNM + $vNM;
	if ($overlaptype eq 'gap') { $totalNM += $overlap; }
	
	# check to see is whole read can be accounted for by alignments to the vector
	my ($isVecRearrange);
	# decide if read is more likely a rearrangement or integration based on edit distance
	# if edit distances are the same, err on the side of caution and assign it as a rearrangement
	$isVecRearrange = isRearrangeOrInt($vCig, $vDir, $vRef, $vPos, $vSec, $vSup, $seq, $thresh, $vNM, $totalNM);

	my ($isHumRearrange);
	$isHumRearrange = isRearrangeOrInt($hCig, $hDir, $hRef, $hPos, $hSec, $hSup, $seq, $thresh, $hNM, $totalNM);
	
	# check to see if location of host alignment is ambiguous: multiple equivalent alignments accounting for host part of read
	my $isHumAmbig;
	if ($hSec eq "NA") { $isHumAmbig = "no";}
	else { $isHumAmbig = isAmbigLoc($hDir, $hCig, $hSec, 'soft', $seq, $tol);}
	
	#check to see if location of viral alignment is ambiguous: multiple equivalent alignments accounting for viral part of read
	my $isVirAmbig;
	if ($vSec eq "NA") { $isVirAmbig = "no";}
	else { $isVirAmbig = isAmbigLoc($vDir, $vCig, $vSec, 'soft', $seq, $tol);}

	### Calculate the start and stop positions of the viral and host sequences relative to the read
	my ($hRStart,$hRStop) = extractSeqCoords($hOri, $hAlig, $overlap, $readlen, $overlaptype);
	my ($vRStart,$vRStop) = extractSeqCoords($vOri, $vAlig, $overlap, $readlen, $overlaptype);
	
	### Collect integration start/stop positions relative to read
	### These are the bases that flank the bond that is broken by insertion of the virus
	### Takes any overlap into account	
	my ($intRStart, $intRStop);

	if    ($hOri eq "+") { ($intRStart,$intRStop) = ($hRStop+1,$vRStart); } 
	elsif ($hOri eq "-") { ($intRStart,$intRStop) = ($vRStop+1,$hRStart); } 
	else 			   { 
		print "Something weird has happened";
		return;
	}
	
	### Assign an orientation - host-virus or virus-host
	my $order;
	# host/virus junction if mapped part is first and read is forward
	# or if mapped part is second and read is reverse
	if    ($hOri eq '+' and $hDir eq 'f') { $order = 'hv'; }
	elsif ($hOri eq '-' and $hDir eq 'r') { $order = 'hv'; }
	
	# vice versa for virus/host
	elsif ($hOri eq '+' and $hDir eq 'r') { $order = 'vh'; }
	elsif ($hOri eq '-' and $hDir eq 'f') { $order = 'vh'; }
	
	
	### Extract viral and host sequence componenets of the read
	### Account for overlap, this won't be included in either the host or viral segments
	
	my $viralSeq = substr($seq, $vRStart-1, $vRStop-$vRStart+1);
	my $hostSeq = substr($seq, $hRStart-1, $hRStop-$hRStart+1);
	my $overlapSeq = substr($seq, $intRStart, $overlap);	
	
	#check that we have extracted the right number of bases (sum of three regions should be whole read)
	if ((length($viralSeq) + length($hostSeq) + length($overlapSeq)) != $readlen) { print "wrong number of bases extracted\n"; }

	### Extract junction coordinates relative to the target sequence
	
	#first get start and stop genomic coordinates of integration
	
	my ($hgStart, $hgStop) = extractCoords($hAlig, $overlap, $overlaptype, $hPos, $order);
	
	# for viral coordinates, we need to pretend the order is opposite (unless virus and host have opposite directions)
	my $viralOrder;
	if ($hDir eq $vDir) {
		if ($order eq 'hv') { $viralOrder = 'vh'; }
		if ($order eq 'vh') { $viralOrder = 'hv'; }
	}
	else {
		$viralOrder = $order;
	}
	my ($vgStart, $vgStop) = extractCoords($vAlig, $overlap, $overlaptype, $vPos, $viralOrder);
	
	my $vIntOri;
	if ($hDir eq 'f') {
		if ($vDir eq 'f') { $vIntOri = '+'; }
		else { $vIntOri = '-'; }
	}
	else {
		if ($vDir eq 'f') { $vIntOri = '-'; }
		else { $vIntOri = '+'; }
	}

	#generate output
	my $outline = join("\t", ($hRef, $hgStart, $hgStop, $vRef, 
				  $vgStart, $vgStop, $overlap, $overlaptype, $order, $vIntOri, 
				  $hostSeq, $viralSeq, $overlapSeq, $hNM, $vNM, $totalNM, $isHumRearrange, 
				  $isVecRearrange, $isHumAmbig, $isVirAmbig, 'chimeric', $hMapq, $vMapq));

}

sub extractOutput {
### Construct the output line
### Output line is a tab seperated string with the following information
### hostChr hostStart hostStop viralChr viralStar viralStop AmbiguousBases OverlapType Orientation(host) hostSeq viralSeq AmbiguousSeq HumanSecAlign ViralSecAlign
	my ($viralData) = shift @_;
	my ($hostData) = shift @_;
	my @intData     = @_;
	
	#viralData and $hostData have the following fields:
	#$refName, $pos, $dir, $cig, $vSec, $vSup, $editDist
	
	#intData has the fields
	#$intRStart, $intRStop, $hRStart, $hRStop, $vRStart, $vRStop, $overlap1, $order, $overlaptype, $isVecRearrange, $isHumRearrange, $isHumAmbig, $isVirAmbig
	
	my @viral = split("\t",$viralData);
	my @host = split("\t",$hostData);

	my $overlap = $intData[6];
	my $overlaptype = $intData[8];
	
	my $vRearrange = $intData[9];
	my $hRearrange = $intData[10];
	
	my $hAmbig = $intData[11];
	my $vAmbig = $intData[12];
	
	my $vNM = $viral[-1];
	my $hNM = $host[-1];
	
	my $totalNM = $hNM + $vNM;
	if ($overlaptype eq 'gap') { $totalNM += $overlap; }

	### Extract junction coordinates relative to the target sequence
	my ($viralStart, $viralStop) = extractCoords($viral[1], $viral[2], $viral[5], $overlap);
	my ($hostStart, $hostStop) = extractCoords($host[1], $host[2], $host[5], $overlap);

	### Extract viral and host sequence componenets of the read
	### Account for overlap, this won't be included in either the host or viral segments
	my $viralSeq = substr($viral[6], $intData[4]-1, ($intData[5]-$intData[4]+1));
	my $hostSeq = substr($host[6], $intData[2]-1, ($intData[3]-$intData[2]+1));

	my $overlapSeq;
	if ($overlap > 0) { # only extract overlap sequence if there is one
		if ($intData[7] eq "hv") { $overlapSeq = substr($viral[6], $intData[3], $overlap); } # order along read is Human -> Viral
		else 			 { $overlapSeq = substr($viral[6], $intData[5], $overlap); } # order along read is Viral -> Human
	}
	else { $overlapSeq = ''; }# leave overlapSeq undefined


	my $outline = join("\t", ($host[0], $hostStart, $hostStop, $viral[0], 
				  $viralStart, $viralStop, $overlap, $overlaptype, $intData[7], 
				  $hostSeq, $viralSeq, $overlapSeq, $hNM, $vNM, $totalNM, $hRearrange, $vRearrange, $hAmbig, $vAmbig, 'chimeric'));

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
		if ($cig =~ /(^\d+)[SH]/)   	  { $ori = "-"; }  # clipped bases first, mapped second
		elsif ($cig =~ /(\d+)[SH]$/) 	  { $ori = "+"; }  # mapped bases first, clipped second
	}
	
	#if read is reverse
	elsif ($dir eq 'r') {
		if ($cig =~ /(^\d+)[SH]/)   	  { $ori = "+"; } 
		elsif ($cig =~ /(\d+)[SH]$/) 	  { $ori = "-"; }

	}
	
	return $ori;
}
