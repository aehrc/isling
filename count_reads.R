#### get reads from folder

library(tidyverse)
library(stringr)
library(writexl)

data_path = "../data/"
files <- list.files(path = data_path, recursive = TRUE, pattern = ".sam")

files

samtools <- "samtools"

#### count mapped and unmapped reads ####

aligns <- data_frame(filename = files) 

aligns <- aligns %>% 
  mutate(mapped = map(filename, ~ system(paste0(samtools, " view -c -F 4 -F 2048 ", getwd(), "/", file.path(data_path, .)), 
                                      intern = TRUE))) %>% 
  unnest()

aligns <- aligns %>% 
  mutate(unmapped = map(filename, ~ system(paste0(samtools, " view -c -f 4 ", getwd(), "/", file.path(data_path, .)), 
                                        intern = TRUE))) %>%
  unnest() 



#add extra columns
aligns <- aligns %>% 
  mutate(mapped = as.numeric(mapped), unmapped = as.numeric(unmapped)) %>% 
  mutate(total = mapped + unmapped) %>% 
  mutate(host_or_virus = ifelse(str_detect(filename, paste0(c("hg38", "mouse"), collapse = "|")), "host", "virus")) %>% 
  mutate(dataset = dirname(dirname(filename))) %>% 
  mutate(sample = str_extract(basename(filename), "^[\\w]+(?=\\.)")) %>%
  mutate(aligned_to = str_extract(basename(filename), "[\\w-]+(?=\\.\\w+\\.sam$)")) %>%
  mutate(align_type = ifelse(str_detect(filename, "bwaPaired"), "paired", "single"))

filtered <- aligns %>% 
  filter((!str_detect(filename, "genomes") & str_detect(dataset, "mouse")))

aligns

#make plot of mapped reads
dsets <- unique(aligns$dataset)
types <- unique(aligns$align_type)
for (j in types) {
	for (i in dsets) {
		aligns %>%
		filter(str_detect(align_type, j)) %>%
  		filter(str_detect(dataset, i)) %>% 
  		gather(key = "type", value = "number", mapped, unmapped) %>%
  		group_by(type, dataset, aligned_to) %>% 
  		ggplot(aes(x = sample, y = number, fill = type)) +
  		geom_bar(stat = "identity", position = "dodge") +
  		theme(axis.text.x=element_text(angle=90, hjust=1)) +
  		facet_wrap(~ aligned_to)
		ggsave(paste0(data_path, i, "/summary/", i, "_alignemnts_", j, ".pdf"))
	}
}

#output xls spreadsheets with numbers of mapped reads
dsets <- unique(aligns$dataset)
types <- unique(aligns$align_type)
for (j in types) {
  for (i in dsets) {
    aligns %>%
      filter(str_detect(align_type, j)) %>%
      filter(str_detect(dataset, i)) %>% 
      gather(key = "type", value = "number", mapped, unmapped) %>%
      group_by(type, dataset, aligned_to) %>% 
      ggplot(aes(x = sample, y = number, fill = type)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme(axis.text.x=element_text(angle=90, hjust=1)) +
      facet_wrap(~ aligned_to)
    ggsave(paste0(data_path, i, "/summary/", i, "_alignments_", j, ".pdf"))
  }
}

for (i in unique(aligns$dataset)) {
  toWrite <- list()
  for (j in unique(aligns$align_type)) {
    toWrite[[j]] <-aligns  %>% 
      filter(str_detect(align_type, j))
  }
  write_xlsx(toWrite, path = paste0(data_path, i, "/summary/", i, "_alignments.xlsx", sep=""))
  
}


#### check fragment sizes (field 9 of sam file) ####


files <- list.files(path = data_path, recursive = TRUE, pattern = "bwaPaired.sam")
#flags: 
# -f 0x40 == only include with flag 0x40 (first in pair)
# -F 0x4 == only exclude with flag 0x4 (read unmapped)
# -F 0x8 == only exclude with flag 0x8 (mate unmapped)
# -f 0x2 == only include with flag 0x2 (mapped in proper pair)

frags <- data_frame(filename = files) %>%
  mutate(proper = map(files, ~ system(paste0(samtools, " view -f 0x40 -f 0x2 ", getwd(), "/", file.path(data_path, .), " | cut -f9"), 
                                      intern = TRUE)))


for (i in files) {
  temp <- frags %>%
    filter(filename == i)
  prop <- tibble(proper_pairs = abs(temp[[1,"proper_pairs"]])) %>%
    filter(proper_pairs > 0 & proper_pairs < thresh)
  ggplot() +
    geom_density(aes(x = prop$proper_pairs)) +
    xlab("Fragment size of proper pairs (bp)")
  ggsave(paste0(data_path, "/", dirname(dirname(i)), "/summary/", basename(i), "_properFrags.pdf"))

}
