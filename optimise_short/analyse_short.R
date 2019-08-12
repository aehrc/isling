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

files <- list.files(data_path, pattern ="insert\\d.short.*txt", recursive=TRUE)

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
  mutate(insPenalty = str_extract(filename, "(?<=insert)\\d+")) %>% 
  mutate(Ithresh = str_detect(filename, "Ithresh"))

#how many conditions had any short integration sites
count(short, total_count > 0 ) %>% filter(`total_count > 0`)
filter(short, total_count > 0)

#calculate the change in number of short integrations for each sample and Ithresh
#ie how does number of short integrations change when we change insPenalty?
short <- short %>%
  group_by(sample, Ithresh) %>%
  mutate(delta_byInsPenalty = total_count - lag(total_count, order_by = insPenalty)) %>%
  mutate(hasChange_byInsPenalty = delta_byInsPenalty !=0) %>%
  ungroup() %>%
  mutate(delta_byIthresh = total_count - lag(total_count, order_by = sample))

#groups where changing insPenalty had an effect on number of short integrations
count(short, delta_byInsPenalty != 0) %>% filter(`delta_byInsPenalty != 0`)

#are there any conditions where there were more short integrations with an Ithresh than without?
#there shouldn't be, because Ithresh should be more stringent

#order in Ithresh column is always True then False for each sample
#so when we calculate the delta above, we alternate between True-False and False-True
#we're looking for positive values for (Ithresh is true) - (Ithresh is false)
#so we want to filter for positive delta_byIthresh and then Ithresh is true
#this table should be empty
filter(short, delta_byIthresh > 0 & Ithresh)


#write table to excel
short %>% 
  select(-short) %>% 
  write_xlsx(path =  paste(out_path, "short_test/count_short.xlsx", sep = ""))


#### plots ####

#plot the number of reads identified by insertion opening penalty
short %>% 
  filter(!Ithresh) %>%
  group_by(sample) %>%
ggplot() +
  geom_line(aes(x = insPenalty, y = delta_byInsPenalty, group = sample)) +
  geom_count(aes(x = insPenalty, y = delta_byInsPenalty)) +
  facet_wrap(vars(dataset)) +
  theme(legend.position='none')
ggsave(filename = paste(out_path, "short_test/count_short_insPenaltydelta.pdf", sep = ""))

short %>% 
  filter(!Ithresh) %>%
  group_by(sample) %>%
ggplot() +
  geom_line(aes(x = insPenalty, y = total_count, group = sample, color=sample)) +
  geom_count(aes(x = insPenalty, y = total_count)) +
  facet_wrap(vars(dataset))  +
  theme(legend.position='none')
ggsave(filename = paste(out_path, "short_test/count_short_insPenalty.pdf", sep = ""), width=20, height=8)
