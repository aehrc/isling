#### combine all integration sites in input files and write to excel spreadsheet ####

#### Packages ####
library(stringr)
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(writexl)
library(tibble)

#### import data ####
args <- commandArgs(trailingOnly=TRUE)

data_files <- args[1:(length(args) - 1)]
out_path <- args[length(args)]
out_path <- file.path(out_path)
out_path <- paste0(out_path, "/")

#import all datasets
df <- tibble::tibble(filename = data_files) %>% # create a data frame holding the file names
  dplyr::mutate(integrations = purrr::map(filename, ~ readr::read_tsv(file.path(.), 
						na = c("", "NA", "?"),
						col_types = cols(Chr = col_character(), NoAmbiguousBases = col_integer(), .default =col_guess())))) %>%
  dplyr::mutate(total_count = flatten_int(map(integrations, ~nrow(.)))) 

#add extra columns with sample name, dataset, host
df <- df %>% 
  dplyr::mutate(dataset = basename(dirname(dirname(filename)))) %>% 
  dplyr::mutate(sample = stringr::str_extract(basename(filename), "^[\\w-_]+(?=\\.)")) %>% 
  dplyr::mutate( host = ifelse(stringr::str_detect(basename(filename), "mouse|mm10|GRCm38"), "mm10", ifelse(stringr::str_detect(basename(filename), "macaque|macaca|macFas5"), "macFas5", "hg38")))

#select all integrations
df <- df %>% 
  dplyr::filter(flatten_lgl(purrr::map(integrations, ~ nrow(.) != 0))) %>% 
  tidyr::unnest(cols = c(integrations))

#write excel spreadsheets for each dataset with one sheet for each sample
for (i in unique(df$dataset)) {
  toWrite <- list()
  data_filt <- df %>% 
    dplyr::filter(dataset == i)
  for (j in unique(data_filt$sample)) {
    toWrite[[j]] <- data_filt  %>% 
      dplyr::filter(sample == j)
  }
  writexl::write_xlsx(toWrite, path = paste(out_path, i, "_annotated.xlsx", sep = ""))
}

# also write an excel spreadsheet for each dataset that doesn't contain the annotations
for (i in unique(df$dataset)) {
  toWrite <- list()
  data_filt <- df %>% 
    dplyr::filter(dataset == i) %>%
    dplyr::select(dataset, sample, host, Chr, IntStart, IntStop, VirusRef, VirusStart, VirusStop, OverlapType, Orientation,
    				HostSeq, ViralSeq, AmbiguousSeq, HostEditDist, ViralEditDist, TotalEditDist, 
    				PossibleHostTranslocation, PossibleVectorRearrangement, HostPossibleAmbiguous, ViralPossibleAmbiguous,
    				Type, ReadID, merged)
  for (j in unique(data_filt$sample)) {
    toWrite[[j]] <- data_filt  %>% 
      dplyr::filter(sample == j)
  }
  writexl::write_xlsx(toWrite, path = paste(out_path, i, ".xlsx", sep = ""))
}

