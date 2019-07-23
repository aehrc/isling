# AAV integration detection pipeline

Inherited from Laurence Wilson.  He wrote it to detect AAV integration sites in data from mice (FRG or wild-type) treated with AAV, from collaborators at Westmead (CMRI).

## Changes from original pipeline

 - Added feature to check for ambiguity in host and virus alignment (ie check for secondary alignments equivalent to the primary one)
 - Added feature to check for possible viral rearrangements or host translocations. If there is one alignment that accounts for the whole read, it's not a possible rearrangement. Otherwise grab all primary, secondary and supplementary alignments and check to see if they can span the whole read collectivley.  If there is a set of alignments that can account for the whole read collectivley, it's a possible vector rearrangement. If there are gaps between adjacent alignments, count the number of bases in those gaps.  If the coverage of the read by alignments is \>95%, it's a possible vector rearrangement, otherwise it's not.
 - Moved most functions from original script to a perl module to facillitate checks for other kinds of integrations: discordant read-pairs, short integrations with (soft-clipped on both ends)
 - Wrote script to check for discordant read-pairs (`discordant.pl`).  Currently this only checks for read-pairs where one read is completely mapped to the host and the other is completely mapped to the virus
 - Implemented check for vector rearrangement as seperate pipeline.  This is relevant for clinical vector where there seems to be a small amount of this going on.  Only need to consider viral alingment here
 - Wrote R scripts for visualization of both integration data and rearrangement data
 - Moved pipeline to snakemake to make it easier to develop
 - In order to save disk space, only consider reads aligned to the virus for alignement to host.  For softclippped reads, extract merged or combined reads that were mapped to the virus, and map to the host.  For discordant read-pairs, extract reads that were mapped to the virus, and reads that didn't map but their mate was mapped.
 - Added de-duplication step prior to alingment (using `clumpify` from `BBMap`).  This can't be done after alignment because different reads might be removed from the host and viral alignments.
