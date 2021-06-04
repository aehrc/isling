# AAV integration detection pipeline

Pipeline to detect viral integration in paired-end reads.

Pipeline requires `snakemake` and either `singularity` (recommended) or `conda` to supply dependencies.  Additionaly, `python` version 3.5 or above and `pandas` are required (these should be automatically installed if installing `snakemake` with `conda`.

Currently have conda environment called `snakemake`, which I'm activating in wrapper script `run_ints.sh`.  This runs the pipeline on the cluster (cluster config `cluster.json`), using conda to fufill dependencies (`envs/*.yml` contains specifications of conda environments).

# Pipeline overview

The pipeline performs several steps in order to identify integration sites.  It takes as input datasets consisting of either fastq files or bam files. It does some pre-processing of the reads (merging overlapping reads, optional) and then aligns them to both a host and a viral sequence.  Reads are first aligned to the viral sequence(s), and then aligned reads are extracted and aligned to the host.  These alignments are used to identify possible viral integrations.

# Running

To run with the test data locally, run:

```
snakemake --configfile test/config/test.yml --cores <cores>
```

## Inputs

A config file, as well as the host and viral references, and reads are required inputs.  Specify all inputs in a config file.

### Config file

A yaml config file is used to specify the data to be processed and parameters for processing.  The config file is organised into datasets, and specifies the parameters to be used for that dataset.  If multiple datasets are present in the config file, they will be processed simultaneously.  Each dataset contains a number of key-value pairs which indicate how the dataset is to be processed.

The first entry in the config file should be the key 'snakefile', with the value being the path to the directory containing the snakefile, for example:

```
snakedir: "/scratch1/sco305/intvi_pipeline"
```

If snakemake is run from the intvi_pipeline directory (which contains the `Snakefile`), this key may be omitted, otherwise it is required.

This should be followed by one block for each dataset to be processed.  All paths should be either absolute or relative to the directory containing the `Snakefile`. An example dataset block from a config file:

```
dataset_name:
  read_folder: "/path/to/read/folder"
  out_dir: "/path/to/output/"
  samples:
  	- "sample1"
  	- "sample2"
  R1_suffix: "_L001_R1.fastq.gz"
  R2_suffix: "_L001_R2.fastq.gz"
  split: 2
  read1-adapt: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
  read2-adapt: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
  mean-frag-len: "estimate"
  dedup: True
  dedup-subs: 2
  merge: True
  trim: True
  host_name: "macFas5"
  host_fasta: "/path/to/refs/macFas5.fa"
  virus_name: "pAAV2-OTC_flop"
  virus_fasta: "/path/to/refs/pAAV2-OTC.fa"
  clip-cutoff: 20
  min-mapq: 10
  cigar-tol: 3
  filter: 
    - "HostEditDist <= 5"
    - "ViralEditDist <= 5"
    - "HostMapQ > 20"
    - "NoAmbiguousBases < 20 or Type == discordant"
  bed-include:
    - "path/to/bed"
  bed-exclude:
    - "path/to/bed"
  min-n-merge: 1
  merge-method: 'common'
```

#### Dataset name

The block should start with a name for the dataset. This can be anything, but each dataset name in the config file must be unique.  This will be used to name output files.

#### Read folder, samples, R1\_suffix, R2\_suffix, bam\_suffix

Specify the path to the folder containing reads for this datset with the `read_folder` key. Reads must be paired-end. Reads can be in either fastq or sam/bam format.

For fastq reads, each R1 and R2 file must have a common suffix. Specify these common suffixes with `R1_suffix` and `R2_suffix`.

For bam files, all reads must be paired and files must have a common suffix, which should be specified with `bam_suffix`.

Specify either the `R1_suffix` and `R2_suffix`, or `bam_suffix`, but not both.

The key `samples` is optional. If it is included, it should be a list of sample names, which together with the suffixes above should describe the filenames of the input reads.  Eg, if the input files are in fastq format, the input `read_folder`, `R1_suffix` and `R2_suffix`, and `samples` list should specify a list of files that can be found at `<read_folder>/<sample><R1_suffix>` and `<read_folder>/<sample><R2_suffix>`.

If the key `samples` is not included, sample names will be inferred from the files present in the `read_folder` having the specified suffixes.  Note that this option requires files to be present on the local filesystem, and is therefore not compatible with cloud execution.

#### Split FASTQ file

The optional `split` feature alllows to split the provided reads in `n` parts to decrease memory usage and to increase amount of possible threads.

#### Adapters and mean fragment length.

Specify the adapters for read 1 and read 2, for trimming and merging, using the `read1-adapt` and `read2-adapt` keys.  Specify either the mean fragment length, or use the string 'estimate' to indicate that the mean fragment length should be estimated from proper pairs in the alignments.

#### Merging, adapter trimming

The keys `merge` and `trim` specify if the reads will be merged (if R1 and R2 are overlapping), or have their adapters trimmed using [SeqPrep](https://github.com/jstjohn/SeqPrep).  The values for these keys must be either `True` or `False`.  If `merge` is `True`, then adapters will be trimmed, regardless of whether `trim` was specified `True` or `False`.

#### Host and virus references

Specify a name and one of either sequence (`fasta` format) or `bwa` index prefix for both the host and virus using the `host_name`, `host_fasta`, `host_prefix`, `virus_name`, `virus_fasta`, `virus_prefix` keys.  The host reference may contian multiple chromosomes, and the virus reference may contain multiple viral sequences.  If both a `fasta` file and `bwa` prefix are specified, only the prefix will be used.

#### De-duplication

The key `dedup` indicates whether duplicate reads should be removed from the reads before merging, trimming and alignment using  [bbmap Clumpify](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/clumpify-guide/).  Sequence-based duplicate removal is necessary since separate alignments are performed to the host and virus, and removing duplicates from these alignments can result in missing reads, since a read pair may appear to be a duplicate in one alignment (and hence would be removed), but may not be a duplicate in the other alignment.  Note that only reads where both members of the pair are duplicates will be removed.  Use the parameter `dedup-subs` to set the number of substitutions allowed for reads to be considered duplicates.  Since this step is memory-intensive, it is not recommended for large datasets.

#### Options for integration detection (`clip-cutoff`, `min-mapq`, `cigar-tol`)

When detecting integrations, only primary alignments with a mapping quality equal to or greater than `min-mapq` will be considered.  This number must be between 0 and 60.  This paramter is optional, and if not set then a default value of 10 will be used.

The parameter `clip-cutoff` is used when detecting chimeric integrations: in the case of a simple junction, this parameter is the minimum nubmer of bases that must come from both host and virus.  That is, for a `clip-cutoff` value of 20, the read must contain at least 20 bases from the host and 20 bases from the virus.  For short integrations (host on either end, virus in the middle), both ends of the read must contain this many bases from the host.  For discordant read pairs, a mapped read may contain at most this number of unmapped bases, and an unmapped read may contain at most this number of mapped bases. This paramter is optional, and if not set then a default value of 20 will be used.

The parameter `cigar-tol` is relevant when there are a small number of non-mapped bases between two mapped regions, or between one mapped region and the end of the read.  If the number of these non-mapped bases is equal to or less than `cigar-tol`, then they will be combined with the nearest mapped region.  For example, for `cigar-tol: 3` `147M3S` would become `150M`, `98M2I50M` would become `150M`, but `145M5S` would remain unchanged.

#### Post-processing

After detection, junction reads/read pairs may be filtered using post-processing.  The following types of post-processing are availble:

1. `filter`: Remove any integrations not meeting user-defined criteria.  Criteria can be based on the following columns:
	- *NoAmbiguousBases* (integer) - the number of bases in a gap or overlap between host and viral alignments
	- *OverlapType* (‘none’, ‘gap’, ‘overlap’, ‘discordant’) - type of overlap between host and viral alignments
	- *Orientation* (‘hv’, ‘vh’) - in the host, is the junction host-virus (+) or virus-host (-)
	- *ViralOrientation* (‘+’, ‘-’) - orientation in which the virus/vector is integrated
	- *HostEditDist* (integer) - edit distance of the alignemnt to the host
	- *ViralEditDist* (integer) - edit distnace of the alignment to the virus/vector
	- *TotalEditDist* (integer) - sum of host and viral edit distances, plus the length of any bases in a gap between host and virus
	- *PossibleHostTranslocation* ('True', 'False') - does the read have two complementary alignments to the host genome?
	- *PossibleVectorRearrangement* ('True', 'False') - does the read have two complementary alignments to the virus/vector
	- *HostAmbiguousLocation* ('True', 'False') - is there a secondary alignment equivalent (same CIGAR) to the primary alignment in the host?
	- *ViralAmbiguousLocation* ('True', 'False') - is there a secondary alignment equivalent (same CIGAR) to the primary alignment in the virus/vector?
	- *Type* (‘chimeric’, ‘discordant’, 'short') - is this read chimeric, a discordant pair, or does it span both junctions (host/virus/host)?
	- *HostMapQ* (integer) - mapping quality of the host alignment
	- *ViralMapQ* (integer) - mapping quality of the virus/vector alignment
	
2. `bed-exclude`: Specify a list of `bed` files in order to exclude any integrations that fall within the regions in those `bed` files
3. `bed-include`: Specify a list of `bed` files in order to only include integrations that fall within the regions in those `bed` files


TODO:
5. `nearest-bed`: Specify a list of `bed` files in order to annotate each integration with the nearest feature in each `bed` file, and the distance between the integration and that feature
6. `nearest-gtf`: Similar to `nearest-bed`, but with a `gtf` file rather than an `bed` file

#### Merging

Integration events with overlapping coordinates and the same orientations are merged.  There are two modes of merging: in `exact` merging, integrations are only merged if their coordinates are exactly the same in both host and virus/vector, and in `common` merging, integrations are merged if they share common coordinates in both host and virus/vector.  The output coordinates in the latter case are the coordinates of the common bases.

Users may be interested only in integrations with a minimum number of supporting reads. The key `merge-n-min` restricts output to integrations with at least this number of integrations.  This file contains the number of reads merged, and the read IDs, so the user can refer back to the summary excel file for more information about each read.

#### Speical dataset options

For convenience, a `global` dataset may additionally be included in the config file.  Options from this dataset will be applied to all other datasets in which the `global` dataset contains an option that is not set in the other dataset.  For keys that are specified in both, the option from the other dataset will be used (ie `global` options don't overwrite values from other datasets).  For example, if the following config file is used:

```
global:
  out_dir: "/path/to/output/"
  R1_suffix: "_L001_R1.fastq.gz"
  R2_suffix: "_L001_R2.fastq.gz"
  read1-adapt: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
  read2-adapt: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
  merge: True
  trim: True
  host_name: "macFas5"
  host_fasta: "/path/to/refs/macFas5.fa"
  virus_name: "pAAV2-OTC_flop"
  virus_fasta: "/path/to/refs/pAAV2-OTC.fa"
  dedup: False
  clip-cutoff: 20
  min-mapq: 10
  cigar-tol: 3
  post:
    - filter
    - dedup
    - mask-exclude:
      - "/path/to/exclude.bed"
    - nearest-gtf:
      - "/path/to/genes.gtf"
 merge-dist: 100
 
dataset-1:
  read_folder: "/path/to/read/folder/1"
  samples:
  	- "sample1"
  	- "sample2"
  
dataset-2:
  read_folder: "/path/to/read/folder/2"
  bam_suffix: ".bam"
```

two datasets (`dataset-1` and `dataset-2`) will be analysed.  `dataset-1` will consist of samples `sample1` and `sample2`, and `dataset-2` will consist of all the `bam` files in the directory `/path/to/read/folder/2`.

### Reads

The reads for each dataset should be placed in the folder specified for each dataset.  Reads must be paired-end.  Note that the pipeline uses SeqPrep to merge and trim reads, and SeqPrep will refuse to process files which contain read-pairs with different read lengths (ie R1 is a different length to R2).  If your reads contain pairs that are different lengths, trim adapters first and then run the pipeline specifing `merge: False` and `trim: False`.

### References

The host and viral references (fasta format) should be specified in the config file for each dataset.  These will be indexed during execution.

## Outputs

A summary of the integrations for each dataset can be found in a folder `summary` in the specified output folder.  For each dataset, the `summary` folder will contain an excel spreadsheet with information about each integration. This file has one tab for each sample in the dataset. Additionally, a `bed` file formatted for use with the [UCSC browswer](https://www.genome.ucsc.edu/) can be found in `summary/ucsc_bed`.

# Types of integrations

## Soft-clipped

Identify soft-clipped reads (`soft.pl`) by identifying reads where one half is mapped and the other soft-clipped in one alignment, and vice versa in the other.  Allow for some gap or overlap between the two alignments.

## Discordant read-pairs

Identify discordant read-pairs (`discordant.pl`) by identifying unmerged pairs where one read is mapped and one unmapped in one alignment, and vice versa for the other alignment.  Allow some soft-clipping in mapped reads and mapping in unmapped reads here - a read is considered mapped if it has a number of soft-clipped bases that is less than a threshould (say 20), and unmapped if it has a number of mapped bases that is less than the same threshould.

## Short integrations

Short integrations are identified on the basis of the following criteria:

1. Host alignment must be matched on both ends, and the length of the match must be more than a cutoff (default 20 base pairs)
2. Host alignment must have an inserted region between the two matched regions
3. Viral alignment must be soft-clipped on both ends with a matched region in the middle.  The number of matched bases must be more than a cutoff (default 20 base pairs)
4. The inserted bases in the host alignment must overlap with the matched bases in the viral alignment (on the read).  This is enforced by checking that the end of the inserted region is more than the start of the matched region, and the end of the matched region is more than the start of the inserted region.

Each read with a short integration has two associated integration sites (one for each end of the integration).

# Identifying ambiguities and possible issues

There are two possible issues that might be occuring.  The first relates to reads that can be accounted for by multiple alignments to the same refrence (possible rearrangement), and the second relates to reads where the primary alignment is equivalent to a secondary alignment, so we can't unambigously assign an integration site (possible location ambiguity).

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

It would be good to distinguish between vector rearrangements and integrations in the regions where the vector and the human genome are homologous (hAAT \[SERPINA1\] promoter, OTC enhancer).  Can do this by looking for variants between vector and human genome in this region: there are a few in the hAAT promoter.  Interrograte this by extracting only the junction fragments from the alignment.  Made a seperate snakefile to do this, in folder `SNPs`.

First pull out reads called as integrations by the pipeline.  Specify regions (usually SERPINA1 promoter and OTC enhancer) in which to call variants and do pileup.  Output vcfs with variants and a pileup.

### Alignment paramters

The pipeline uses bwa-mem for alignment.  The settings are more or less default, except that the mismatch penalty has been lowered from 4 to 2.  I think this might be leading to fewer intergration sites, because sub-optimal alignments are scored higher than better ones that might lead to missing some integration sites.

One example of this is read M01996:226:000000000-CFBJG:1:2107:24838:18761 in sample 22 of the london macaque data.  This read is not indicated to be an integration site. CIGAR in chr1 is 97M159H (supplementary alignment, reverse), POS is 120875219, main alignment is Chr7 CIGAR 101S31M8I48M2I66M (reverse).  Viral alignment is forward 156M100S.  Two primary alignments don't give integration because clipped part is on the same side of the read, but secondary host alignment with viral primary could be integration. Host coordinates of matched region would be chr1:12075122-120875219 - this is 297 bp from Sunando's site.

In general, optimising the alingment might be difficult because a) we don't really know what a true positive in this datset is and b) optimising parameters might lead to overfitting.
