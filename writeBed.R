#### this script takes the output from the viral integration pipeline, and outputs a bed file for display in the UCSC browser
#### it filters the data according to the criteria below:
#### - exclude vector rearrangments
#### - exclude host translocations
#### - exclude any integrations with a gap or overlap of more than 20 bp (unless they are discordant)
#### - exclude any integrations with a host or virus edit distance of more than 5 bp

#### load libraries ####
library(stringr)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(readr)

#### path to the data ###
# pass in files to process as command line arguments
# within this folder, each dataset is stored in a seperate folder
# within the dataset folders, the integrations for each sample are contained in a folder called ints, 
# 


args <- commandArgs(trailingOnly=TRUE)

data_files <- args[1:(length(args) - 1)]

out_path <- args[length(args)]
out_path <- file.path(out_path)
out_path <- paste0(out_path, "/")

dir.create(file.path(out_path), showWarnings = FALSE)



#import all datasets
df <- tibble::tibble(filename = data_files) %>% # create a data frame holding the file names
  dplyr::mutate(integrations = purrr::map(filename, ~ readr::read_tsv(file.path(.), 
						na = c("", "NA", "?"),
						col_types = cols(Chr = col_character(), NoAmbiguousBases = col_integer(), .default =col_guess())))) %>%
  dplyr::mutate(total_count = purrr::flatten_int(purrr::map(integrations, ~nrow(.)))) 


#add extra columns with sample name, dataset, host
df <- df %>% 
  dplyr::mutate(dataset = dirname(dirname(filename))) %>% 
  dplyr::mutate(sample = stringr::str_extract(basename(filename), "^[\\w-_]+(?=\\.)")) %>% 
  dplyr::mutate( host = dplyr::case_when(
				stringr::str_detect(basename(filename), "mouse|mm10|GRCm38") ~ "mm10",
				stringr::str_detect(basename(filename), "macaque|macaca|macFas5") ~ "macFas5",
				stringr::str_detect(basename(filename), "hg19") ~ "hg19",
				stringr::str_detect(basename(filename), "hg38|GRCh38|human") ~ "hg38",
				TRUE ~ "???"
))


#unnest
df <- df %>% 
  dplyr::filter(purrr::flatten_lgl(map(integrations, ~ nrow(.) != 0))) %>% 
  tidyr::unnest(integrations)



# write bed files for importing into UCSC browser
# minimally required fields are Chr, Start and Stop
# chr field must begin with "chr"

#note that need to add "chr" to chromosome numbers
#chrom chromStart chromEnd

chroms <- paste0("chr", c(1:22, "X", "Y", "M")) %>% paste(collapse = "|")


for (i in unique(df$dataset)) {
  toWrite <- list()
  data_filt <- df %>% 
    dplyr::filter(dataset == i)
  for (j in unique(data_filt$sample)) {
    df %>%
	dplyr::filter(sample == j) %>%
        dplyr::mutate(Chr = ifelse(Chr == "MT", "M", Chr)) %>%
	dplyr::mutate(Chr = ifelse(sum(stringr::str_detect(Chr, "chr")), Chr, paste0("chr", Chr))) %>%
	dplyr::select(Chr, IntStart, IntStop, ReadID) %>%
	dplyr::filter(str_detect(Chr, chroms)) %>%
	dplyr::rename(`#chrom` = Chr) %>%
	dplyr::rename(`ChromStart` = IntStart) %>%
	dplyr::rename(`ChromEnd` = IntStop) %>%
	dplyr::rename(name = ReadID) %>%
	readr::write_tsv(path = paste0(out_path, basename(i), ".", j, ".post.bed"), col_names = TRUE)
  }
}



