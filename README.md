# AAV integration detection pipeline

Inherited from Laurence Wilson.  He wrote it to detect AAV integration sites in data from mice (FRG or wild-type) treated with AAV, from collaborators at Westmead (CMRI).

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

## To do

 - Write a third script, `short.pl` to look for short insertions.  These would appear as a mapped portion on both sides in the human alignment with an insertion in the middle, and in the viral alingment they would appear as soft-clipped on both sides with a mapped region in the middle.
 - Check errors in `discordant.pl`
 - Fix bug in `softClip.pl` where sometimes can't get host/viral sequences from the read
