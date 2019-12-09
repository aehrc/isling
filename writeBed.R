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
# data are contained in a folder ../out
# within this folder, each dataset is stored in a seperate folder
# within the dataset folders, the integrations for each sample are contained in a folder called ints, 
# 


data_path = "../out/"

out_path = "../out/summary/ucsc_bed/"
dir.create(file.path(out_path), showWarnings = FALSE)

#get integration files
files <- list.files(data_path, pattern =".integrations.txt", recursive=TRUE)

#import all datasets
df <- tibble(filename = files) %>% # create a data frame holding the file names
  mutate(integrations = map(filename, ~ read_tsv(file.path(data_path, .), 
						na = c("", "NA", "?"),
						col_types = cols(Chr = col_character(), NoAmbiguousBases = col_integer(), .default =col_guess())))) %>%
  mutate(total_count = flatten_int(map(integrations, ~nrow(.)))) 

#add extra columns with sample name, dataset, host
df <- df %>% 
  mutate(dataset = dirname(dirname(filename))) %>% 
  mutate(sample = str_extract(basename(filename), "^[\\w]+(?=\\.)")) %>% 
  mutate( host = ifelse(str_detect(basename(filename), "mouse|mm10|GRCm38"), "mm10", ifelse(str_detect(basename(filename), "macaque|macaca|macFas5"), "macFas5", "hg38")))

#unnest
df <- df %>% 
  filter(flatten_lgl(map(integrations, ~ nrow(.) != 0))) %>% 
  unnest(integrations)


#get 'clean' integrations':

#change names of datasets
df <- df %>% 
  mutate(dataset = case_when(
    str_detect(dataset, "FRG_OTC") ~ "FRG1",
    str_detect(dataset, "AGRF_CAGRF21377_CNW6N") ~ "FRG2",
    str_detect(dataset, "macaque") ~ "macaque",
    TRUE ~ dataset
  ))

#defne regions to annotate for each host
hosts <- c("macFas5", "macFas5", "hg38", "hg38")
feats <- c("OTC", "SERPINA1", "OTC", "SERPINA1")
chroms <- c("X", "7", "X", "14")
starts <- c(36515003, 159438228, 38343260, 94388561)
stops <- c(36515195, 159440447, 38343435, 94390786)

#label the data
df <- df %>% 
  mutate(label = 
           case_when(
             host == hosts[1] & Chr == chroms[1] & ((IntStart > starts[1] & IntStart < stops[1]) | (IntStop > starts[1] & IntStop < stops[1])) ~ feats[1],
             host == hosts[2] & Chr == chroms[2] & ((IntStart > starts[2] & IntStart < stops[2]) | (IntStop > starts[2] & IntStop < stops[2])) ~ feats[2],
             host == hosts[3] & Chr == chroms[3] & ((IntStart > starts[3] & IntStart < stops[3]) | (IntStop > starts[3] & IntStop < stops[3])) ~ feats[3],
             host == hosts[4] & Chr == chroms[4] & ((IntStart > starts[4] & IntStart < stops[4]) | (IntStop > starts[4] & IntStop < stops[4])) ~ feats[4],
             TRUE ~ "other"
           )
  )

#filter integrations that don't meet our criteria
clean <- df %>% 
  filter(HostEditDist <= 5) %>% 
  filter(ViralEditDist <= 5) %>% 
  filter(NoAmbiguousBases <= 20 | OverlapType == 'discordant') %>% 
  filter(label == 'other') %>% 
  filter(PossibleVectorRearrangement == 'no') %>% 
  filter(PossibleHostTranslocation == "no") 


# write bed files for importing into UCSC browser
# minimally required fields are Chr, Start and Stop
# chr field must begin with "chr"

#note that need to add "chr" to chromosome numbers
#chrom chromStart chromEnd

chroms <- paste0("chr", c(1:22, "X", "Y", "M")) %>% paste(collapse = "|")


for (i in unique(clean$dataset)) {
  toWrite <- list()
  data_filt <- clean %>% 
    filter(dataset == i)
  for (j in unique(data_filt$sample)) {
    clean %>%
	filter(sample == j) %>%
	mutate(Chr = paste0("chr", Chr)) %>%
	select(Chr, IntStart, IntStop, ReadID) %>%
	filter(str_detect(Chr, chroms)) %>%
	rename(`#chrom` = Chr) %>%
	rename(`ChromStart` = IntStart) %>%
	rename(`ChromEnd` = IntStop) %>%
	rename(name = ReadID) %>%
	write_tsv(path = paste0(out_path, i, ".", j, ".clean.bed"), col_names = TRUE)
  }
}



