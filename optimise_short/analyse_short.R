#### Packages ####
library(stringr)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(readr)
library(writexl)

#### import data ####
data_path = "../../out/"
out_path = "../../out/summary/"

files <- list.files(data_path, pattern ="insert\\d.short.txt", recursive=TRUE)

files

#import all datasets
short <- tibble(filename = files) %>% # create a data frame holding the file names
  mutate(short = map(filename, ~ read_tsv(file.path(data_path, .), 
                                                 na = c("", "NA", "?")))) %>%
  mutate(total_count = flatten_int(map(short, ~nrow(.))))

#add extra columns with metadata extracted from filename
short <- short %>% 
  mutate(dataset = dirname(dirname(filename))) %>% 
  mutate(sample = str_extract(basename(filename), "^[\\w]+(?=\\.)")) %>% 
  mutate( host = ifelse(str_detect(basename(filename), "mouse|mm10"), "mm10", "hg38")) %>% 
  mutate(insPenalty = str_extract(filename, "(?<=insert)\\d+"))

#write table to excel
short %>% 
  select(-short) %>% 
  write_xlsx(path =  paste(out_path, "count_short.xlsx", sep = ""))


#### plots ####

#plot the number of reads identified by insertion opening penalty
short %>% 
ggplot() +
  geom_line(aes(x = insPenalty, y = total_count, group=sample)) +
  geom_count(aes(x = insPenalty, y = total_count))
ggsave(filename = paste(out_path, "count_short_insPenalty.pdf", sep = ""))
