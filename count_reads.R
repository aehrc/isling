#### get reads from folder

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(writexl)

data_path = "../out/"

#### count mapped and unmapped reads ####

aligns <- read_tsv("../out/summary/count_mapped.txt")

#add extra columns
aligns <- aligns %>% 
  mutate(total = mapped + unmapped) %>%
  mutate(sample_short = str_extract(sample, "^[[:alnum:]\\-_]+(?=\\.)")) %>% 
  mutate(aligned_to = str_match(sample, "^[[:alnum:]\\-_]+\\.([[:alnum:]\\-_]+?)\\.")[,2]) %>% 
  mutate(align_type =  str_match(sample, "^[[:alnum:]\\-_]+?\\.[[:alnum:]\\-_]+?\\.([[:alpha:]\\-_]+?)\\.")[,2]) %>% 
  select(-total) #total column is incorrect

# add column for if alignment is viral or not == doesn't contain "Mappedreads"
#also for junction reads == contain "juncts"
aligns <- aligns %>% 
  mutate(viralAlign = !str_detect(sample, "Mappedreads")) %>% 
  mutate(junctAlign = str_detect(sample, "juncts"))

write_xlsx(aligns, path = paste0(data_path, "summary/count_mapped.xlsx"))

#calculate total reads for each sample - 
#this is the number of mapped and unmapped reads for the viral alignment
#because when doing host alignment only used reads that were mapped to the virus

juncts <- aligns %>% 
  filter(junctAlign)

aligns <- aligns %>% 
  filter(!junctAlign)

#this only works for the virus alignments
aligns <- aligns %>% 
  mutate(total = case_when(
    viralAlign ~ mapped + unmapped,
    TRUE ~ as.integer(0)
  ))

#calculate the fraction of the total reads aligned to the virus/vector
aligns %>% 
  mutate(frac_aligned = mapped / total) %>% 
  filter(viralAlign) %>% 
  filter(align_type == "bwa") %>% 
  select(dataset, sample, aligned_to, mapped, unmapped, frac_aligned) %>% 
  write_xlsx(path = paste0(data_path, "summary/frac-mapped-to-virus.xlsx"))

#make plot of mapped reads
dsets <- unique(aligns$dataset)
types <- unique(aligns$align_type)
for (j in types) {
	for (i in dsets) {
		aligns %>%
		filter(align_type == j) %>%
  		filter(dataset == i) %>% 
  		gather(key = "type", value = "number", mapped, unmapped) %>%
  		group_by(type, dataset, aligned_to) %>% 
  		ggplot(aes(x = sample_short, y = number, fill = type)) +
  		geom_bar(stat = "identity", position = "dodge") +
  		theme(axis.text.x=element_text(angle=90, hjust=1)) +
  		facet_wrap(~ aligned_to)
		ggsave(paste0(data_path, "summary/aligns_", i, "_", j, ".pdf", sep=""))
	}
}

#output xls spreadsheets with numbers of mapped reads
for (i in unique(aligns$dataset)) {
  toWrite <- list()
  for (j in unique(aligns$align_type)) {
    toWrite[[j]] <-aligns  %>% 
      filter(align_type == j) %>% 
      filter(dataset == i)
  }
  write_xlsx(toWrite, path = paste0(data_path, "summary/aligns_", i, ".xlsx", sep=""))
  
}


