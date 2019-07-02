mo#### get reads from folder

library(tidyverse)

align_path = "data/2019-06-17_alignments/"
files <- dir(align_path, pattern =".sam")

samtools <- "/Users/sco305/Documents/Projects/Integration/expt2_pipeline-tweaks/tools/samtools-1.9/samtools"

aligns <- data_frame(filename = files) %>%
  mutate(mapped = map(files, ~ system(paste0(samtools, " view -c -F 4 -F 2048 ", getwd(), "/", file.path(align_path, .)), 
                                      intern = TRUE, 
                                      timeout = 30))) %>% 
  unnest() %>% 
  mutate(unmapped = map(files, ~ system(paste0(samtools, " view -c -f 4 ", getwd(), "/", file.path(align_path, .)), 
                                        intern = TRUE, 
                                        timeout = 30))) %>%
  unnest() 

#add extra columns
aligns <- aligns %>% 
  mutate(mapped = as.numeric(mapped), unmapped = as.numeric(unmapped)) %>% 
  mutate(total = mapped + unmapped) %>% 
  mutate(type = ifelse(str_detect(filename, paste0(c("hg38", "mouse"), collapse = "|")), "host", "virus")) %>% 
  mutate( dataset = ifelse(str_detect(filename, "FRG"), "FRG", 
                           ifelse(str_detect(filename, "combinedReads"), "patient", 
                                  ifelse(str_detect(filename, "m3"), "mouse", "second_patient")))) %>% 
  mutate(sample = ifelse(str_detect(dataset, "mouse"), substr(filename, 0, 5), substr(filename, 0, 3))) 

filtered <- aligns %>% 
  filter((!str_detect(filename, "genomes") & str_detect(dataset, "mouse")))

aligns

#make plot of mapped reads
colns = c("FRG", "patient", "mouse")
for (i in colns) {
aligns %>%
    filter(str_detect(dataset, i)) %>% 
  gather(key = "type", value = "number", mapped, unmapped) %>%
  group_by(type, dataset) %>% 
  ggplot(aes(x = filename, y = number, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x=element_text(angle=90, hjust=1)) 
ggsave(paste0("out/reads_", i, ".pdf"))
}


