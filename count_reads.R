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
  mutate(sample_short = str_match(aligns$sample, "^([[:alnum:]-_]+)\\.[[:alnum:]-_]+\\.[:alpha:]+\\.[[:alnum:]-_]+\\.bam")[,2]) %>% 
  mutate(aligned_to = str_match(aligns$sample, "^[[:alnum:]-_]+\\.([[:alnum:]-_]+)\\.[:alpha:]+\\.[[:alnum:]-_]+\\.bam")[,2]) %>% 
  mutate(align_type = str_match(aligns$sample, "^[[:alnum:]-_]+\\.[[:alnum:]-_]+\\.([:alpha:]+)\\.[[:alnum:]-_]+\\.bam")[,2])

write_xlsx(aligns, path = paste0(data_path, "summary/count_mapped.xlsx"))

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


