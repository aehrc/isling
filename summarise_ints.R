#### Plot integration site FRG206 ####

#### Packages ####
library(stringr)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(readr)
library(writexl)
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

#check rate of false negatives for vector rearrangement:
#look for events on chromosome X or 14 in mouse data aligned against human
#that aren't marked as possible vector rearrangmenets

mouse_align_hg38 <- data %>% 
  filter(dataset == "mouse_with_OTC") %>% 
  filter(host == "hg38") %>% 
  mutate(on_14_or_x = str_detect(Chr, "14|X"))

#get number of 'false negatives' in possible vector rearrangement
mouse_align_hg38 %>% 
  filter(on_14_or_x) %>% 
  group_by(PossibleVectorRearrangement) %>% 
  summarise(vecRearrange = n())

#write this data to excel spreadsheet
toWrite <- list()
for (i in unique(mouse_align_hg38$sample)) {
  toWrite[[i]] <-mouse_align_hg38 %>% 
    filter(str_detect(sample, i))
}
write_xlsx(toWrite, path =  paste(out_path, "mouse_with_OTC_hg38align.xlsx", sep = ""))

#remove mouse data aligned against human from the rest of analysis
data <- data %>% 
filter(!(dataset == "mouse_with_OTC" & host == "hg38"))

data %>%
  group_by(dataset, PossibleVectorRearrangement, host) %>% 
  summarise( count = n())


#write excel spreadsheets for each dataset with one sheet for each sample
for (i in unique(data$dataset)) {
  toWrite <- list()
  data_filt <- data %>% 
    filter(dataset == i)
  for (j in unique(data_filt$sample)) {
    toWrite[[j]] <- data_filt  %>% 
      filter(sample == j)
  }
  write_xlsx(toWrite, path = paste(out_path, i, ".xlsx", sep = ""))
}

#get summary of number of sites per dataset
data %>% 
  group_by(dataset, filename) %>% 
  count()  %>% 
  write_xlsx(path = paste0(out_path, "num_sites.xlsx"))

#loop over columns to make pie charts of all columns
colns <- c("OverlapType", "PossibleVectorRearrangement", "PossibleHostTranslocation", "HostPossibleAmbiguous", "ViralPossibleAmbiguous", "has_issue", "has_rearrangement", "has_locAmbig")
for (i in colns) {
  data %>%
    mutate(dataset = as.factor(dataset)) %>% #need to make this a factor for complete to work
    count(!!ensym(i), dataset) %>% # get counts by coln[i] and dataset
    complete(!!ensym(i), dataset, fill = list(n = 0)) %>% #fill in zero rows
    group_by(dataset) %>% #need this to get frequency of counts by dataset
      mutate( freq = n/sum(n) ) %>% #get frequency
      ggplot(aes(x="", y=freq, fill=!!ensym(i))) +
      geom_bar(stat = "identity", width=1) +
      geom_label(aes(label = paste0(n, "\n", round(freq*100), "%")), 
                 position = position_fill(vjust = 0.5),
                 show.legend = FALSE) +
      coord_polar("y", start=0) +
      theme_void() +
      theme(legend.position = "bottom") +
      facet_wrap(vars(dataset))
  ggsave(paste0(out_path, "pie_",gsub(" ", "", i), ".pdf", sep = ""))
}


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

