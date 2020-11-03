### post-process integrations from the perl scripts
### specify the file to be processed and the types of processing to do on it

## type of of post-processing:
## types will be applied in the order in which they appear below, regardless of which types are specified

### dedup: remove integrations with identical read sequences - keep only one representative integration
### filter: filter integrations according to a set of criteria (ie gap length < 20bp, NM <5, etc etc)
### mask-exclude: remove integrations that fall within regions specified in a bed file (specify the bed file)
### mask-include: only include integrations that fall within regions specified in a bed file
### nearest-gtf: annotate the nearest genomic feature from a gtf file
### nearest-bed: annotate the nearest genomic feature from a bed file

## after post-processing, files will be saved in their original location
## with ".postprocessed" added before the file extension

#### packages ####
library(dplyr)
library(readr)
cat("\n")

#### parse command-line arguments ####
# first command line argument should be path to file to postprocess
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1)
{
  stop("provide path to file(s to process")
}

if (length(args) < 2)
{
  cat("no types specified!\n")
}

# types of post-processing available
types <- c("filter", "dedup", "mask-exclude", "mask-include", "nearest-gtf", "nearest-bed", "RNA-seq-gtf")

# check which types of post-processing to do
# check if we want to do de-duplication
dedup <- FALSE
if ("dedup" %in% args)
{
  dedup <- TRUE
}
# check if we want to do filtering
filt <- FALSE
if ("filter" %in% args)
{
  filt <- TRUE
}
# check if we want to exclude based on a bed file
mask_exclude <- ""
if ("mask-exclude" %in% args)
{
  # get bed files to use for masking
  if (which(args == "mask-exclude") + 1 > length(args))
  {
    stop("specify bed files containing regions to use for masking")
  }
  
  mask_exclude <- args[which(args == "mask-exclude") + 1]
  
  # check that we didn't get any other types (indicating a missing bed file)
  if (sum(types %in% mask_exclude) > 0)
  {
    stop("specify bed files containing regions to use for masking")
  }
  
  # check that file extension of specified bed files is ".bed"
  for (i in mask_exclude)
  {
    extension <- tools::file_ext(i)
    if (extension[length(extension)] != "bed")
    {
      stop("specify bed files containing regions to use for masking")
    }
  }
  
  # check that files exist
  if (sum(file.exists(mask_exclude)) != length(mask_exclude))
  {
    stop("(some) bed file(s) do not exist")
  }

}
# check if we want to only include integration in a bed file
mask_include <- ""
if ("mask-include" %in% args)
{
  # get bed files to use for masking
  if (which(args == "mask-include") + 1 > length(args))
  {
    stop("specify bed files containing regions to use for masking")
  }
  
  mask_include <- args[which(args == "mask-include") + 1]
  
  # check that we didn't get any other types (indicating a missing bed file)
  if (sum(types %in% mask_include) > 0)
  {
    stop("specify bed files containing regions to use for masking")
  }
  
  # check that file extension of specified bed files is ".bed"
  for (i in mask_include)
  {
    extension <- tools::file_ext(i)
    if (extension[length(extension)] != "bed")
    {
      stop("specify bed files containing regions to use for masking")
    }
  }
  
  # check that files exist
  if (sum(file.exists(mask_include)) != length(mask_include))
  {
    stop("(some) bed file(s) do not exist")
  }
  
}

# check if we want to annotate based on a gtf file
nearest_gtf <- ""
if ("nearest-gtf" %in% args)
{
  # get bed files to use for masking
  if (which(args == "nearest-gtf") + 1 > length(args))
  {
    stop("specify gtf files for annotation")
  }
  
  nearest_gtf <- args[which(args == "nearest-gtf") + 1]
  
  # check that we didn't get any other types (indicating a missing bed file)
  if (sum(types %in% nearest_gtf) > 0)
  {
    stop("specify gtf files for annotation")
  }
  
  # check that file extension of specified gtf files is ".gtf"
  for (i in nearest_gtf)
  {
    extension <- tools::file_ext(i)
    if (extension[length(extension)] != "gtf")
    {
      stop("specify gtf files for annotation")
    }
  }
  
  # check that files exist
  if (sum(file.exists(nearest_gtf)) != length(nearest_gtf))
  {
    stop("(some) gtf file(s) do not exist")
  }
}

# check if we want to annotate based on a bed file
nearest_bed <- ""
if ("nearest-bed" %in% args)
{
  # get bed files to use for masking
  if (which(args == "nearest-bed") + 1 > length(args))
  {
    stop("specify bed files for annotation")
  }
  
  nearest_bed <- args[which(args == "nearest-bed") + 1]
  
  # check that we didn't get any other types (indicating a missing bed file)
  if (sum(types %in% nearest_bed) > 0)
  {
    stop("specify bed files for annotation")
  }
  
  # check that file extension of specified bed files is ".bed"
  for (i in nearest_bed)
  {
    extension <- tools::file_ext(i)
    if (extension[length(extension)] != "bed" & extension[length(extension)] != "tsv")
    {
      stop("specify bed/tsv files for annotation")
    }
  }
  
  # check that files exist
  if (sum(file.exists(nearest_bed)) != length(nearest_bed))
  {
    stop("(some) bed file(s) do not exist")
  }
}

# check if we want to use RNA-seq
nearest_gtf_RNA <- ""
if ("RNA-seq-gtf" %in% args)
{
  # get gtf files to use for nearest
  if (which(args == "RNA-seq-gtf") + 1 > length(args))
  {
    stop("specify gtf files for annotation")
  }
  
  # need to provide the GTF file as the first argument after the RNA-seq-gtf argument
  nearest_gtf_RNA <- args[which(args == "RNA-seq-gtf") + 1]
  print(nearest_gtf_RNA)
  
  # get the name of the column to get from the RNA-seq tsv files
  RNA_seq_col <- args[which(args == "RNA-seq-gtf") + 2]
  
  # get the tsv file containing RNA-seq counts
  RNA_seq_tsv <- args[which(args == "RNA-seq-gtf") + 3]

  #check that all these vectors are the same length
  if ((length(nearest_gtf_RNA) != length(RNA_seq_col)) | (length(RNA_seq_col) != length(RNA_seq_tsv))) {
    stop("specify a gtf file, column and tsv file for each RNA-seq annotation")
  }
  
  # check that we didn't get any other types (indicating a missing bed file)
  if (sum(types %in% nearest_gtf) > 0)
  {
    stop("specify gtf files for annotation")
  }
  
  # check that file extension of specified gtf files is ".gtf"
  for (i in nearest_gtf_RNA)
  {
    extension <- tools::file_ext(i)
    if (extension[length(extension)] != "gtf")
    {
      stop("specify gtf files for annotation")
    }
  }
  
  # check that files exist
  if (sum(file.exists(nearest_gtf_RNA)) != length(nearest_gtf_RNA))
  {
    stop("(some) gtf file(s) do not exist")
  }
  
  # check that files exist
  if (sum(file.exists(RNA_seq_tsv)) != length(RNA_seq_tsv))
  {
    stop("(some) RNA-seq txt file(s) do not exist")
  }
}


#### data import ####
# columns are dependent on the output of the perl scripts that identify integrations
int_cols <- readr::cols(
  Chr = col_character(),
  IntStart = col_integer(),
  IntStop = col_integer(),
  VirusRef = col_character(),
  VirusStart = col_integer(),
  VirusStop = col_integer(),
  NoAmbiguousBases = col_integer(),
  OverlapType = col_character(),
  Orientation = col_character(),
  HostSeq = col_character(),
  ViralSeq = col_character(),
  AmbiguousSeq = col_character(),
  HostEditDist = col_integer(),
  ViralEditDist = col_integer(),
  TotalEditDist = col_integer(),
  PossibleHostTranslocation = col_character(),
  PossibleVectorRearrangement = col_character(),
  HostPossibleAmbiguous = col_character(),
  ViralPossibleAmbiguous = col_character(),
  Type = col_character(),
  ReadID = col_character(),
  merged = col_character()
)

# read integration file
ints <- readr::read_tsv(args[1],
                 col_types = int_cols,
                 na = c("?", "", "NA"))

cat(nrow(ints), "integrations imported\n")

#### de-duplicating ####

if (dedup & nrow(ints) != 0)
{
  source("post/dedup.R")
}

#### filtering ####

if (filt & nrow(ints) != 0 & nrow(ints) != 0)
{
  source("post/filter.R")
}

#### masking ####

# do exluding based on intersect with bed file
if (mask_exclude[1] != "" & nrow(ints) != 0)
{
  source("post/mask-exclude.R")
}

# do excluding based on lack of intersect with bed file
if (mask_include[1] != "" & nrow(ints) != 0)
{
  source("post/mask-include.R")
}

#### annotating ####

# annotate nearest feature from gtf file(s)
if (nearest_gtf[1] != "" & nrow(ints) != 0)
{
  source("post/nearest-gtf.R")
}

#annotate nearest feature from bed file(s)
if (nearest_bed[1] != "" & nrow(ints) != 0)
{
  source("post/nearest-bed.R")
}

#annotate nearest feature from bed file(s)

if (nearest_gtf_RNA[1] != "" & nrow(ints) != 0)
{
  source("post/RNA-seq-gtf-tsv.R")
}

#### save output postprocessed file ####

base_name <- tools::file_path_sans_ext(args[1])
ext <- tools::file_ext(args[1])

readr::write_tsv(ints, path = paste0(base_name, ".post.", ext))


#### session info ####
cat("\n")
#sessionInfo()
