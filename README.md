# AAV integration detection pipeline

Inherited from Laurence Wilson.  He wrote it to detect AAV integration sites (chimeric reads only) in data from mice (FRG or wild-type) treated with AAV, from collaborators at Westmead (CMRI).

Pipeline requires conda and snakemake.  Currently have conda environment called `snakemake`, which I'm activating in wrapper script `run_ints.sh`.  This runs the pipeline on the cluster (cluster config `cluster.json`), using conda to fufill dependencies (`envs/*.yml` contains specifications of conda environments).

## Changes from original pipeline

 - Added feature to check for ambiguity in host and virus alignment (ie check for secondary alignments equivalent to the primary one)
 - Added feature to check for possible viral rearrangements or host translocations. If there is one alignment that accounts for the whole read, it's not a possible rearrangement. Otherwise grab all primary, secondary and supplementary alignments and check to see if they can span the whole read collectivley.  If there is a set of alignments that can account for the whole read collectivley, it's a possible vector rearrangement. If there are gaps between adjacent alignments, count the number of bases in those gaps.  If the coverage of the read by alignments is \>95%, it's a possible vector rearrangement, otherwise it's not.
 - Moved most functions from original script to a perl module to facillitate checks for other kinds of integrations: discordant read-pairs, short integrations with (soft-clipped on both ends)
 - Wrote script to check for discordant read-pairs (`discordant.pl`).  This checks for read pairs that have one read mapped and one unmapped in the host, and the converse in a virus.  It makes use of a paired alignment, so the pipeline now does both a paired (with all the pairs) and a single (with the merged and the unmerged reads) alignment.  The perl script then parses both the viral and host paired alingments and looks for pairs meeting the the criteria.  It allows for some soft-clipping that's less than a cutoff (20 bp, same as the one used for `softclip.pl`).
 - Moved pipeline to snakemake to make it easier to develop.
 - Implemented check for vector rearrangement as seperate pipeline.  This is relevant for clinical vector where there seems to be a small amount of this going on.  Only need to consider viral alingment here.  Snakefile is `rearrange.sf`, wrapper script is `run_rearrange.sh`.
 - Wrote R scripts for visualization of both integration data `summarise_ints.R` and rearrangement data `count_rearrange.R`.
 - In order to save disk space, only consider reads aligned to the virus for alignement to host.  For soft-clippped reads, extract merged or combined reads that were mapped to the virus, and map to the host.  For discordant read-pairs, extract reads that were mapped to the virus, and reads that didn't map but their mate was mapped.
 - Added de-duplication step prior to alingment (using `clumpify` from `BBMap`).  This can't be done after alignment because different reads might be removed from the host and viral alignments.
 - Added script `short.pl` to look for short insertions.
	- Wrote script to first identify reads that look like they might be short insertions: clipped on both sides in viral alignment (more than cutoff bases), mapped on both ends in human alignment (more than cutoff bases) with insertion in the middle
	- Used snakefile and scripts in `optimise_short` to try to optimise the alignment to identify more short insertions.  First tried to vary the penatly for a new insertion between 0 and the default (6).  See section optimise_short below.

## Short insertions

Short insertions are identified on the basis of the following criteria:
1. Host alignment must be matched on both ends, and the length of the match must be more than a cutoff (default 20 base pairs)
2. Host alignment must have an inserted region between the two matched regions
3. Viral alignment must be soft-clipped on both ends with a matched region in the middle.  The number of matched bases must be more than a cutoff (default 20 base pairs)
4. The inserted bases in the host alignment must overlap with the matched bases in the viral alignment (on the read).  This is enforced by checking that the end of the inserted region is more than the start of the matched region, and the end of the matched region is more than the start of the inserted region.

Each read with a short insertion has two associated integration sites (one for each end of the insertion).  

### Optimize\_short

This folder contains scripts for optimization of identifying short insertions.  The idea was to try different alignment parameters in order to find the parameters that give the most insertions.  Tried this on a number of different datasets: most CMRI datasets, a subest (100 samples) of human/hbv dataset (PRJNA298941), as well as mouse/AAV dataset (PRJNA485509).

Varied penalty for creating an insertion between 0 and 6.  In general, most samples had no integration sites, some had one or two.  The most were found in the human/hbv dataset, although inspection of these indicated that within each sample, mostly it's just one site with lots of supporting reads.  Judging by number of short integrations, it's best to have a penalty of 0.  Howver, this does sometimes result in CIGAR containing a few smaller inserted regions seperated by multiple small matched regions, rather than just one larger insertion.  If going this route, might need to combine adjacent M/I CIGAR operations in decection script.

### Location of short insertions

One issue that came up is that for short insertions we can see both sides of the inserted viral region.  So we'd expect both to be bookended in the human genome.  However, in the first attempt at characterising these events, quirks in the alignment meant that they were not.  This comes about because of CIGARS like 58M5I7M14I7M12I16M5I112M in the host.  There are multiple small matched ares breaking up the inserted region.  For the purposes of short integrations, I treat the whole middle region, from the first to the last inserted bases, as one big insertion.  The overlaps/gaps and ambiguous bases are calculated based on this big insertion.  

However, when calculating the genomic coordinates, I use the CIGAR and the 1-based leftmost mapping position (from the SAM file) to calculate the coordinates for each CIGAR operation, and then pick the first and last matched region for the host coordinates.  Therefore, any matched regions in the middle will throw off the genomic coordinates.

There are two ways to get around this:
1. Adjust alignment - make the gap creation penalty more than the gap extension penalty.  Hopefully this results in one long insertion, rather than one broken up by multiple matches.  This could be fiddly
2. Combine all inserted regions in CIGAR within script.  This might be easier because it doesn't involve mucking around with alignment.  Try this first.



## To do

 - Write a third script, `short.pl` to look for short insertions.  These would appear as a mapped portion on both sides in the human alignment with an insertion in the middle, and in the viral alingment they would appear as soft-clipped on both sides with a mapped region in the middle.
 - Check errors in `discordant.pl`
 - Fix bug in `softClip.pl` where sometimes can't get host/viral sequences from the read
 - Get genetic elements in which integrations occur
 - Improved visualisations
