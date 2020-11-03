# this script annotates integration sites based on their proximity to genomic features
# Suzanne Scott, 11/03/2020
# suzanne.scott@csiro.au

# genes/exons/transcripts/etc should be provided as a gtf file containing only the features of interest
# other genomic features can be provided as bed files.  Calls to `bedtools nearest` is used to annotate the
# nearest feature and the distance to it.

# pass in a tsv file containing integrations.  
# this must contain the columns Chr and IntStart, specifying the chromosome and 
# position on that chromosome of each integration.  
# It must also have a  
# Other metadata may also be present

#### libraries ####
library(tidyverse)

#### parse command line arguments ####
# first command line argument should be path to integrations.txt file that 
# comes out of the perl scripts identifying integrations
args = commandArgs(trailingOnly=TRUE)

ints <- readr::read_tsv(file = args[1])

#https://stackoverflow.com/questions/14525580/how-to-pass-command-line-arguments-when-calling-source-on-an-r-file-within-ano
