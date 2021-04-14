#### combine all integration sites in input files and write to excel spreadsheet ####

# usage: Rscirpt summarise_ints.R host virus merged1 merged2 ... mergedn outdir

#### Packages ####
library(stringr)
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(writexl)
library(tibble)
library(tools)
library(glue)

#### import data ####
args <- commandArgs(trailingOnly=TRUE)
host <- args[1]
virus <- args[2]
data_files <- args[3:(length(args) - 1)]
out_path <- args[length(args)]
out_path <- file.path(out_path)
out_path <- paste0(out_path, "/")

# function to get samples from data_files
getSamples <- function(data_files) {
	data_files <- basename(data_files)
	data_files <- data_files %>%
									file_path_sans_ext() %>%
									file_path_sans_ext() %>%
									file_path_sans_ext() %>%
									file_path_sans_ext() %>%
									file_path_sans_ext() 
	return(data_files)
}


#import all datasets
df <- tibble::tibble(filename = data_files) %>% # create a data frame holding the file names
  dplyr::mutate(integrations = purrr::map(filename, ~ readr::read_tsv(file.path(.), 
						na = c("", "NA", "?"),
						col_types = cols(Chr = col_character(), .default =col_guess())))) 

#add extra columns with sample name, dataset, host
df <- df %>% 
  dplyr::mutate(dataset = basename(dirname(dirname(filename)))) %>% 
  dplyr::mutate(sample = stringr::str_match(basename(filename), glue("(.+)\\.{host}\\.{virus}\\.integrations\\.post\\.merged\\.txt"))[,2])

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
  fn <- paste(out_path, i, "_annotated.xlsx", sep = "")
  print(glue("saving to {fn}"))
  writexl::write_xlsx(toWrite, path = fn)
}

# also write an excel spreadsheet for each dataset that doesn't contain the annotations
if (nrow(df) == 0) {
	dset <- basename(dirname(dirname(data_files[1])))
	samples <- getSamples(data_files)
	toWrite <- list()
	for (i in samples) {
		toWrite[[i]] <- tibble()
	}
	fn1 <- paste(out_path, dset, ".xlsx", sep = "")
	fn2 <- paste(out_path, dset, "_annotated.xlsx", sep = "")
	print(glue("saving to {fn1}"))
	print(glue("saving to {fn2}"))
	writexl::write_xlsx(toWrite, path = fn1)
	writexl::write_xlsx(toWrite, path = fn2)
} else {
	for (i in unique(df$dataset)) {
	  toWrite <- list()
	  data_filt <- df %>% 
		    dplyr::filter(dataset == i) %>%
    dplyr::select(sample, Chr, IntStart, IntStop, Virus, VirusStart, VirusStop, nChimeric, nDiscordant, SiteID, ReadIDs)
	  for (j in unique(data_filt$sample)) {
	    toWrite[[j]] <- data_filt  %>% 
	      dplyr::filter(sample == j) %>% 
	      dplyr::select(-sample)
	  }
	  fn <- paste(out_path, i, ".xlsx", sep = "")
	  print(glue("saving to {fn}"))
	  writexl::write_xlsx(toWrite, path = fn)
	}
}

sessionInfo()

