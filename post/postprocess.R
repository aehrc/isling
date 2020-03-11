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

#### packages ####
library(tidyverse)
cat("\n")

#### parse command-line arguments ####
# first command line argument should be path to file to postprocess
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 1)
{
  stop("provide path to file to process")
}
if (length(args) < 2)
{
  cat("no types specified!\n")
}

# types of post-processing available
types <- c("filter", "mask-exclude", "mask-include", "nearest-gtf", "nearest-bed")

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
    extension <- stringr::str_split(i, "\\.")[[1]]
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
    extension <- stringr::str_split(i, "\\.")[[1]]
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

#### data import ###
# columns are dependent on the output of the perl scripts that identify integrations
int_cols <- readr::cols(
  Chr = col_character(),
  IntStart = col_double(),
  IntStop = col_double(),
  VirusRef = col_character(),
  VirusStart = col_double(),
  VirusStop = col_double(),
  NoAmbiguousBases = col_double(),
  OverlapType = col_character(),
  Orientation = col_character(),
  HostSeq = col_character(),
  ViralSeq = col_character(),
  AmbiguousSeq = col_character(),
  HostEditDist = col_double(),
  ViralEditDist = col_double(),
  TotalEditDist = col_double(),
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

if (dedup)
{
  source("dedup.R")
}

#### filtering ####

if (filt)
{
  source("filter.R")
}

#### masking ####
 
if (mask_exclude[1] != "")
{
  source("mask-exclude.R")
}

if (mask_include[1] != "")
{
  source("mask-include.R")
}

#### annotating ####







#### session info ####
cat("\n")
sessionInfo()
