#### this script takes the output from the viral integration pipeline, and outputs a bed file for display in the UCSC browser

# usage: Rscript writeBed.R host virus [<bed1> <bed2> ... <bedn>] <out_path>

#### load libraries ####
library(stringr)
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(glue)

#### path to the data ###
# pass in files to process as command line arguments
# within this folder, each dataset is stored in a seperate folder
# within the dataset folders, the integrations for each sample are contained in a folder called ints, 
# 

args <- commandArgs(trailingOnly=TRUE)
host <- args[1]
virus <- args[2]
data_files <- args[3:(length(args) - 1)]

out_path <- args[length(args)]
out_path <- file.path(out_path)
out_path <- paste0(out_path, "/")

dir.create(file.path(out_path), showWarnings = FALSE)



#import all datasets
int_cols <- readr::cols(
  Chr = col_character(),
  IntStart = col_integer(),
  IntStop = col_integer(),
  JunctionType = col_character(),
  Virus = col_character(),
  VirusStart = col_integer(),
  VirusStop = col_integer(),
  VirusOri = col_character(),
  nChimeric = col_integer(),
  nDiscordant = col_integer(),
  SiteID = col_integer(),
  ReadIDs = col_character(),
  .default = col_guess()
)

df <- tibble::tibble(filename = data_files) %>% # create a data frame holding the file names
  dplyr::mutate(integrations = purrr::map(filename, ~ readr::read_tsv(file.path(.), 
						na = c("", "NA", "?"),
						col_types = int_cols)))  

#add extra columns with sample name, dataset, host
df <- df %>% 
  dplyr::mutate(dataset = dirname(dirname(filename))) %>% 
  dplyr::mutate(sample = stringr::str_match(basename(df$filename), glue("(.+)\\.{host}\\.{virus}.integrations.post.merged.txt"))[,2])



#unnest
df <- df %>% 
  dplyr::filter(purrr::flatten_lgl(map(integrations, ~ nrow(.) != 0))) %>% 
  tidyr::unnest(integrations)



# write bed files for importing into UCSC browser
# minimally required fields are Chr, Start and Stop
# chr field must begin with "chr"

#note that need to add "chr" to chromosome numbers
#chrom chromStart chromEnd

if (nrow(df) == 0) {
	exp <- basename(out_path)
	file.create(paste0(out_path, "empty.bed"))
} else {
	for (i in unique(df$dataset)) {
	  toWrite <- list()
	  data_filt <- df %>% 
		dplyr::filter(dataset == i)
	  for (j in unique(data_filt$sample)) {
	  df %>%
		dplyr::filter(sample == j) %>%
		dplyr::mutate(Chr = dplyr::case_when(
			!stringr::str_detect(Chr, "chr") ~ paste0("chr", Chr),
			Chr == "chrMT" ~ "chrM",
			TRUE ~ Chr
			)) %>%
	  dplyr::mutate(Chr = ifelse(Chr == "chrMT", "chrM", Chr)) %>%
		dplyr::select(Chr, IntStart, IntStop, ReadIDs) %>%
		dplyr::filter(stringr::str_detect(Chr, 'chr\\d+|chrX|chrY|chrM')) %>%
		dplyr::rename(`#chrom` = Chr) %>%
		dplyr::rename(`ChromStart` = IntStart) %>%
		dplyr::rename(`ChromEnd` = IntStop) %>%
		dplyr::rename(name = ReadIDs) %>%
		readr::write_tsv(path = paste0(out_path, basename(i), ".", j, ".post.bed"), col_names = TRUE)
	  }
   }
}


sessionInfo()
