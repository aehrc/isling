#### Plot integration site FRG206 ####

#### Packages ####
library(stringr)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(readr)
library(circlize)
 
#### import data ####
data_path = "../out/"
out_path = "../out/summary/"


files <- list.files(data_path, pattern =".integrations.txt", recursive=TRUE)

#remove empty files
isnotempty <- map(files, ~file.info(paste(data_path, ., sep=''))$size != 0) %>% unlist()
files <- files[isnotempty]

#import all datasets
data <- data_frame(filename = files) %>% # create a data frame holding the file names
  mutate(integrations = map(filename, ~ read_tsv(file.path(data_path, .), 
						na = c("", "NA", "?"), 
						col_types = cols(Chr = col_character(), NoAmbiguousBases = col_integer(), .default =col_guess())))) %>% 
  unnest()


#add extra columns with sample name, dataset, host
data <- data %>% 
  mutate(dataset = dirname(dirname(filename))) %>% 
  mutate(sample = str_extract(basename(filename), "^[\\w]+(?=\\.)")) %>% 
  mutate( host = ifelse(str_detect(basename(filename), "mouse|mm10"), "mm10", "hg38"))


#add summary columns based on ambiguites/issues
data <- data %>% 
  mutate(has_issue = ifelse((PossibleVectorRearrangement == "yes" | 
                               PossibleHostTranslocation == "yes" |
                               HostPossibleAmbiguous == "yes" |
                               ViralPossibleAmbiguous == "yes"), 'yes', 'no')) %>% 
  mutate(has_rearrangement = ifelse(PossibleVectorRearrangement == "yes" | 
                                      PossibleHostTranslocation == "yes", "yes", "no")) %>% 
  mutate(has_locAmbig = ifelse( HostPossibleAmbiguous == "yes" |
                                  ViralPossibleAmbiguous == "yes", 'yes', 'no'))



#discard PCR duplicates: filter for unique read sequences (merged for softClip, R1+R2 for discordant)
data <- data %>% 
  distinct(merged, .keep_all = TRUE)


#remove mouse data aligned against human from the rest of analysis
data <- data %>% 
filter(!(dataset == "mouse_with_OTC" & host == "hg38"))


#make circos plots for all datasets
datasets <- unique(data$dataset)
for (i in datasets) {
  
  #get host for this dataset
  host <- data %>% 
    filter(str_detect(dataset, i)) %>% 
    select(host) %>% 
    unique()
  
  #get chromosomes for host
  chroms <- if (str_detect(host, "hg38")) {paste0("chr", c(1:22, "X", "Y"))} else {paste0("chr", c(1:20, "X", "Y"))}
  
  #make bed file of data with no issues
  bed_noissues <- data %>% 
    filter(str_detect(dataset, i)) %>% 
    filter(has_issue == "no") %>% 
    select(Chr:IntStop) %>%
    mutate(Chr = paste0("chr", Chr)) %>% 
    filter(str_detect(Chr, paste(chroms, collapse = "|"))) %>% 
    as.data.frame()
  
  #make bed file of data with issues
  bed_issues <- data %>% 
    filter(str_detect(dataset, i)) %>% 
    filter(has_issue == "yes") %>% 
    select(Chr:IntStop) %>%
    mutate(Chr = paste0("chr", Chr)) %>% 
    filter(str_detect(Chr, paste(chroms, collapse = "|"))) %>% 
    as.data.frame()
  
  #make bed file of all data
  bed <- data %>% 
    filter(str_detect(dataset, i)) %>% 
    select(Chr:IntStop) %>%
    mutate(Chr = paste0("chr", Chr)) %>% 
    filter(str_detect(Chr, paste(chroms, collapse = "|"))) %>% 
    as.data.frame()
  
  #make plot of all data
  pdf(paste0(out_path, "hostCircos_", i, "_density_all.pdf", sep = ""))
  if (str_detect(host, "hg38")) {circos.initializeWithIdeogram(species = host, 
                                                               chromosome.index = chroms)} else {circos.initializeWithIdeogram(species = host)}
  circos.genomicRainfall(bed, track.height = 0.2)
  circos.genomicDensity(bed, track.height = 0.1)
  dev.off()
  circos.clear()
  
  #make plot of data with no issues
  bed_list2 = list(bed_issues, bed_noissues)
  pdf(paste0(out_path, "hostCircos_", i, "_density_issues.pdf", sep = ""))
  if (str_detect(host, "hg38")) {circos.initializeWithIdeogram(species = host, 
                                                               chromosome.index = chroms)}
  else {circos.initializeWithIdeogram(species = host)}
  circos.genomicRainfall(bed_list2, col = c("#00BFC4", "#F8766D"), track.height = 0.2)
  circos.genomicDensity(bed_issues, col = c("#00BFC4"), track.height = 0.1)
  circos.genomicDensity(bed_noissues, col = c("#F8766D"), track.height = 0.1)
  dev.off()
  circos.clear()

}









