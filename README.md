# AAV integration detection pipeline

Inherited from Laurence Wilson.  He wrote it to detect AAV integration sites (chimeric reads only) in data from mice (FRG or wild-type) treated with AAV, from collaborators at Westmead (CMRI).

Pipeline requires conda and snakemake.  Currently have conda environment called `snakemake`, which I'm activating in wrapper script `run_ints.sh`.  This runs the pipeline on the cluster (cluster config `cluster.json`), using conda to fufill dependencies (`envs/*.yml` contains specifications of conda environments).


# Pipeline overview

The pipeline performs several steps in order to identify integration sites.  It takes as input datasets consisting of either fastq files or bam files. It does some pre-processing of the reads (de-duplication and merging overlapping reads, both optional) and then aligns them to both a host and a viral sequence.  These alignments are used to identify possible viral integrations.

## Inputs

### Config file

A config file is used to specify the data to be processed and parameters for processing.  The config file is organised into datasets, and specifies the paramters to be used for that dataset.  If multiple datasets are present in the config file, they will be processed simultaneously.  

An example dataset from a config file:
```
dataset_name:
  host: "GRCh38"
  virus: "AAV_genomes"
  merge: "True"
  dedup: "True"
  read1-adapt: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
  read2-adapt: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
  R1_suffix: "_L001_R1.fastq.gz"
  R2_suffix: "_L001_R2.fastq.gz"
  post:
    - filter
    - dedup
```

### Reads

Data should be paired-end.  Each dataset is processed with the same parameters (host, virus, postprocessing), and should be kept together in the same folder.



# Types of integrations

## Soft-clipped

Identify soft-clipped reads (`soft.pl`) by identifying reads where one half is mapped and the other soft-clipped in one alignment, and vice versa in the other.  Allow for some gap or overlap between the two alignments.
	
## Discordant read-pairs

Identify discordant read-pairs (`discordant.pl`) by identifying unmerged pairs where one read is mapped and one unmapped in one alignment, and vice versa for the other alignment.  Allow some soft-clipping in mapped reads and mapping in unmapped reads here - a read is considered mapped if it has a number of soft-clipped bases that is less than a threshould (say 20), and unmapped if it has a number of mapped bases that is less than the same threshould.

## Short insertions

Short insertions are identified on the basis of the following criteria:
1. Host alignment must be matched on both ends, and the length of the match must be more than a cutoff (default 20 base pairs)
2. Host alignment must have an inserted region between the two matched regions
3. Viral alignment must be soft-clipped on both ends with a matched region in the middle.  The number of matched bases must be more than a cutoff (default 20 base pairs)
4. The inserted bases in the host alignment must overlap with the matched bases in the viral alignment (on the read).  This is enforced by checking that the end of the inserted region is more than the start of the matched region, and the end of the matched region is more than the start of the inserted region.

Each read with a short insertion has two associated integration sites (one for each end of the insertion).  

### Optimize\_short

This folder contains scripts for optimization of identifying short insertions.  The idea was to try different alignment parameters in order to find the parameters that give the most insertions.  Tried this on a number of different datasets: most CMRI datasets, a subest (100 samples) of human/hbv dataset (PRJNA298941), as well as mouse/AAV dataset (PRJNA485509).

Varied penalty for creating an insertion between 0 and 6.  In general, most samples had no integration sites, some had one or two.  The most were found in the human/hbv dataset, although inspection of these indicated that within each sample, mostly it's just one site with lots of supporting reads.  Judging by number of short integrations, it's best to have a penalty of 0.  However, this does sometimes result in CIGAR containing a few smaller inserted regions seperated by multiple small matched regions, rather than just one larger insertion.  

Decided to combine small insertion/deletions/other CIGAR elements together with neighbouring elements.  Needed different way of doing this than LW came up with in `processCIGAR`.  In the context of soft-clipped integration sites, he picked out the matched part of the read and then combined everything else and called it soft-clipped.  This won't work if I'm looking for short integrations. Instead, combine any CIGAR elements that are less than a certain length (defined as `$tol` with a default of 3).  Combine only if the operations either side are the same, and if the operation to be combined consumes the query then add the bases to the combined operation, otherwise ignore them. 

For example:

`50S20M2I20M10S`: the operation `2I` has neighbours `20M` and `20M`.  Its length is less than `$tol`, so it will be combined. The `I` CIGAR operation consumes the query so its bases are added to the combined operation.  The final CIGAR is `50S42M10S`

`50S20M2D20M10S`: the operation `2D` has neighbours `20M` and `20M`.  Its length is less than `$tol`, so it will be combined. The `D` CIGAR operation does not consume the query so its bases are not added to the combined operation.  The final CIGAR is `50S40M10S`

`50S20M2I10S`: the operation `2I` is short, but the operations on either side (`20M`, `10S`) are different.  The CIGAR will be unchanged

`50S20M5I20M10S`: the operation `5I` has neighbours `20M` and `20M`, but its length is more than `$tol`, so it will not be combined. The CIGAR will be unchanged.

### Location of short insertions

One issue that came up is that for short insertions we can see both sides of the inserted viral region.  So we'd expect both to be bookended in the human genome.  However, in the first attempt at characterising these events, quirks in the alignment meant that they were not.  This comes about because of CIGARS like 58M5I7M14I7M12I16M5I112M in the host.  There are multiple small matched ares breaking up the inserted region.  For the purposes of short integrations, I treat the whole middle region, from the first to the last inserted bases, as one big insertion.  The overlaps/gaps and ambiguous bases are calculated based on this big insertion.

However, when calculating the genomic coordinates, I use the CIGAR and the 1-based leftmost mapping position (from the SAM file) to calculate the coordinates for each CIGAR operation, and then pick the first and last matched region for the host coordinates.  Therefore, any matched regions in the middle will throw off the genomic coordinates.

There are two ways to get around this:
1. Adjust alignment - make the gap creation penalty more than the gap extension penalty.  Hopefully this results in one long insertion, rather than one broken up by multiple matches.  This could be fiddly
2. Combine all inserted regions in CIGAR within script.  This might be easier because it doesn't involve mucking around with alignment.  Try this first.

# Identifying ambiguities and possible issues

Claus identified two possible issues that might be occuring.  The first relates to reads that can be accounted for by multiple alignments to the same refrence (possible rearrangement), and the second relates to reads where the primary alignment is equivalent to a secondary alignment, so we can't unambigously assign an integration site (possible location ambiguity).

## Location ambiguity

To identify reads where a location of the integration can't be unambigouosly defined, look for alignments where the CIGAR for the primary alignment is the same as the CIGAR for any of the secondary alignemnts.  Note that this approach disregards mismatches, which are appear as mapped in the CIGAR.

## Possible rearrangements

 For the clincal OTC vector datasets, we see a lot of intergration in regions of the host genome that are contained in the vector.  The clinical vector contains the hAAT (SERPINA1) promoter and the OTC enhancer, and there were a lot of apparent integrations in these regions in the host genome.  However, a lot of these reads which were indicated to be integrations could also be completely accounted for by alignments to the vector genome, indicating that they could also be vector rearrangments.  Translocations are also possible, where the read is accounted for by multiple alingments to the host genome.

At first, identified possible rearrangements (either vector rearrangements or host translocations) as reads where multiple alignments to one genome span most of the read.  A first attempt at this check was to see if 95% of bases in read were accounted for by any number of alignments from one genome, that read is flagged as a possible rearrangement

However, this led to some rearrangments being missed - an analysis of the variants in these regions, particularly in the FRG and macaque datasets, indicated that there are very few, if any, true integrations in these regions.  I therefore changed this check to consider two explanations for every read that might be an integration:
 - The read is an integration, composed of viral and host portions
 - The read is a vector rearragement (or a host translocation), composed of multiple alignments to only one reference genome
For each possible explanation, I calculate an 'edit distance', which is defined as the sum of the edit distances for the individual alingments, and the number of bases in any gaps between those alignments.  The edit distance for each of these explanations is then compared, and the explanation with the smaller edit distance is output.  In the case of a tie, vector rearrangement is output to err on the side of caution.

### Variants

It would be good to distinguish between vector rearrangements and insertions in the regions where the vector and the human genome are homologous (hAAT \[SERPINA1\] promoter, OTC enhancer).  Can do this by looking for variants between vector and human genome in this region: there are a few in the hAAT promoter.  Interrograte this by extracting only the junction fragments from the alignment.  Made a seperate snakefile to do this, in folder `SNPs`.

First pull out reads called as integrations by the pipeline.  Specify regions (usually SERPINA1 promoter and OTC enhancer) in which to call variants and do pileup.  Output vcfs with variants and a pileup.


### Alignment paramters

The pipeline uses bwa-mem for alignment.  The settings are more or less default, except that the mismatch penalty has been lowered from 4 to 2.  I think this might be leading to fewer intergration sites, because sub-optimal alignments are scored higher than better ones that might lead to missing some integration sites.

One example of this is read M01996:226:000000000-CFBJG:1:2107:24838:18761 in sample 22 of the london macaque data.  This read is not indicated to be an integration site. CIGAR in chr1 is 97M159H (supplementary alignment, reverse), POS is 120875219, main alignment is Chr7 CIGAR 101S31M8I48M2I66M (reverse).  Viral alignment is forward 156M100S.  Two primary alignments don't give integration because clipped part is on the same side of the read, but secondary host alignment with viral primary could be integration. Host coordinates of matched region would be chr1:12075122-120875219 - this is 297 bp from Sunando's site.

In general, optimising the alingment might be difficult because a) we don't really know what a true positive in this datset is and b) optimising parameters might lead to overfitting.



## To do


 - Fix bug in `softClip.pl` where sometimes can't get host/viral sequences from the read
 - Improve check for vector rearrangement - still a lot of integrations in SERPINA1/OTC that aren't marked as possible vector rearrangements
 - Annotate integrations - genes/promoters/etc that they appear in
 - Simulate integrations and reads to 
 - Try out structural variation tools as another approach
 - Improved visualisations
 
# Changes from original pipeline

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
	- Used snakefile and scripts in `optimise_short` to try to optimise the alignment to identify more short insertions.  First tried to vary the penatly for a new insertion between 0 and the default (6).  See section optimise\_short.

