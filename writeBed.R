#### this script takes the output from the viral integration pipeline, and outputs a bed file for display in the UCSC browser
#### it filters the data according to the criteria below:
#### - exclude vector rearrangments
#### - exclude host translocations
#### - exclude any integrations with a gap or overlap of more than 20 bp (unless they are discordant)
#### - exclude any integrations with a host or virus edit distance of more than 5 bp

#### load libraries ####
library(tidyverse)

#### path to the data ###
# data are contained in a folder ../out
# within this folder, each dataset is stored in a seperate folder
# within the dataset folders, the integrations for each sample are contained in a folder called ints, 
# 
data_path = "../out/"



