#### Packages ####
library(stringr)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(readr)
library(writexl)

#import data

data_path = "../../out/"
out_path = "../../out/summary/rearrange/"
refiles <- list.files(path = data_path, recursive = TRUE, pattern = "rearrange.txt")

#remove empty files
isnotempty <- map(refiles, ~file.info(paste(data_path, ., sep=''))$size != 0) %>% unlist()
refiles <- refiles[isnotempty]

colns <- c("ReadID", "possibleVecRearrange", "totalGapBP", "no_gaps", "sups", "seq")

#import all datasets
rearrange <- data_frame(filename = refiles) 

rearrange <- rearrange %>% # create a data frame holding the file names
  mutate(rearrange = map(filename, ~ read_tsv(file.path(data_path, .), 
                                                 na = c("", "NA", "?"), 
                                                 col_names = colns,
                                              skip = 1))) %>% 
  unnest()



#add extra columns with sample name, dataset, host
rearrange <- rearrange %>% 
  mutate(dataset = dirname(dirname(filename))) %>% 
  mutate(sample = str_extract(basename(filename), "^[\\w]+(?=\\.)")) %>% 
  mutate( host = ifelse(str_detect(basename(filename), "mouse|mm10"), "mm10", "hg38"))




#make summary table of number of vector rearrangements
suumary <- rearrange %>% 
  group_by(possibleVecRearrange, dataset) %>% 
  count() %>%
  write_xlsx(path = paste0(out_path, "num_rearrange.xlsx", sep = ""))

#make pie chart of possible vector rearrange
colns <- c("no_gaps","possibleVecRearrange")
for (i in colns) {
  rearrange %>% 
  mutate(dataset = as.factor(dataset)) %>%
  count(!!ensym(i), dataset) %>% 
  complete(!!ensym(i), dataset, fill = list(n = 0))%>% 
    group_by(dataset) %>% 
    mutate( freq = n/sum(n) ) %>% 
    ggplot(aes(x="", y=freq, fill=!!ensym(i))) +
    geom_bar(stat = "identity", width=1)+
    geom_label(aes(label = paste0(n, "\n", round(freq*100), "%")), 
               position = position_fill(vjust = 0.5),
               show.legend = FALSE) +
    coord_polar("y", start=0) +
    theme_void() +
    theme(legend.position = "bottom") + 
    facet_wrap(vars(dataset))
  ggsave(paste0(out_path, "rearrange_", i, ".pdf", sep = ""))

} #end of loop over columns
  

#clear arrange
remove(rearrange)


#make location histogram for viral rearrangements
#load bed file
bedfiles <- list.files(path = data_path, recursive = TRUE, pattern = "rearrange.bed")

#remove empty files
isnotempty <- map(bedfiles, ~file.info(paste(data_path, ., sep=''))$size != 0) %>% unlist()
bedfiles <- bedfiles[isnotempty]



#import all datasets
colns <- c("virus", "start", "stop", "rearrange", "ID", "seq")
bed <- data_frame(filename = bedfiles) %>% # create a data frame holding the file names
  mutate(integrations = map(filename, ~ read_tsv(file.path(data_path, .), 
                                                 na = c("", "NA", "?"),
                                                 col_names = colns))) %>% 
  unnest()

#calculate midpoint of bed

bed <- bed %>% 
  mutate(dataset = dirname(dirname(filename))) %>% 
  mutate(sample = str_extract(basename(filename), "^[\\w]+(?=\\.)")) %>% 
  mutate( host = ifelse(str_detect(basename(filename), "mouse|mm10"), "mm10", "hg38"))


bed <-bed %>% 
  mutate(mid = (start + stop)/2)


summary <- bed %>% 
  filter(!str_detect(filename, "AGRF")) %>% 
  filter(str_detect(sample, "m")) %>% 
  group_by(sample, rearrange) %>% 
  count()

bed %>% 
  filter(rearrange == "yes") %>% 
  filter(str_detect(filename, "pAAV")) %>% 
  ggplot(aes(x = mid, color = sample))+ 
  geom_density() +
  facet_grid(rows = vars(dataset))
ggsave(file.path(out_path, "rearrange_density.pdf"), width = 10, height = 5)
  
