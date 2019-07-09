#### Packages ####
library(tidyverse)

#import data

data_path = "../data/"
out_path = "../out/"
files <- list.files(path = data_path, recursive = TRUE, pattern = "rearrange.txt")

#remove empty files
isnotempty <- map(files, ~file.info(paste(data_path, ., sep=''))$size != 0) %>% unlist()
files <- files[isnotempty]


#import all datasets
data <- data_frame(filename = files) %>% # create a data frame holding the file names
  mutate(integrations = map(filename, ~ read_tsv(file.path(data_path, .), 
                                                 na = c("", "NA", "?"), 
                                                 col_types = cols(Chr = col_character(), .default =col_guess())))) %>% 
  unnest()


#add extra columns with sample name, dataset, host
data <- data %>% 
  mutate(dataset = dirname(dirname(filename))) %>% 
  mutate(sample = str_extract(basename(filename), "^[\\w]+(?=\\.)")) %>% 
  mutate( host = ifelse(str_detect(basename(filename), "mouse|mm10"), "mm10", "hg38"))




#make summary table of number of vector rearrangements
data %>% 
  group_by(possibleVecRearrange) %>% 
  count()

#make pie chart of possible vector rearrange
colns <- c("possibleVecRearrange", "Gaps")
for (i in colns) {
  count(data, !!ensym(i)) %>% 
    mutate( freq = n/sum(n) ) %>% 
    ggplot(aes(x="", y=freq, fill=!!ensym(i))) +
    geom_bar(stat = "identity", width=1)+
    geom_label(aes(label = paste0(n, "\n", round(freq*100), "%")), 
               position = position_fill(vjust = 0.5),
               show.legend = FALSE) +
    coord_polar("y", start=0) +
    theme_void() +
    theme(legend.position = "bottom")
  ggsave(paste0(out_path, "rearrange_",gsub(" ", "", i), ".pdf", sep = ""))

} #end of loop over columns
  


#make location histogram for viral rearrangements
#load bed file
bedfiles <- list.files(path = data_path, recursive = TRUE, pattern = "rearrange.bed")

#remove empty files
isnotempty <- map(bedfiles, ~file.info(paste(data_path, ., sep=''))$size != 0) %>% unlist()
bedfiles <- bedfiles[isnotempty]


#import all datasets
data <- data_frame(filename = bedfiles) %>% # create a data frame holding the file names
  mutate(integrations = map(filename, ~ read_tsv(file.path(data_path, .), 
                                                 na = c("", "NA", "?"), 
                                                 col_types = cols(Chr = col_character(), .default =col_guess())))) %>% 
  unnest()


  