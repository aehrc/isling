library(tidyverse)

align_path = "test_data/"
files <- dir(align_path, pattern ="Paired.sam")
thresh <- 1000

samtools <- "/Users/sco305/Documents/Projects/Integration/expt2_pipeline-tweaks/tools/samtools-1.9/samtools"

frags <- tibble(filename = files) %>% 
  mutate(all_pairs = map(files, ~as.numeric(system(paste0(samtools, " view -f 0x40 -F 0x4 -F 0x8 ", getwd(), "/", file.path(align_path, .), " | cut -f9"), intern = TRUE)))) %>% 
  mutate(proper_pairs = map(files, ~as.numeric(system(paste0(samtools, " view -f 0x40 -f 0x2 ", getwd(), "/", file.path(align_path, .), " | cut -f9"), intern = TRUE)))) 

frags %>% 
  mutate(all_pairs = as.vector(all_pairs))

for (i in files) {
  temp <- frags %>% 
    filter(filename == i) 
  prop <- tibble(proper_pairs = abs(temp[[1,3]])) %>% 
    filter(proper_pairs > 0 & proper_pairs < thresh)
  ggplot() +
    geom_density(prop, aes(x = proper_pairs)) +
    xlab("Fragment size of proper pairs (bp)")
}

#flags: 
# -f 0x40 == only include with flag 0x40 (first in pair)
# -F 0x4 == only exclude with flag 0x4 (read unmapped)
# -F 0x8 == only exclude with flag 0x8 (mate unmapped)
# -f 0x2 == only include with flag 0x2 (mapped in proper pair)

all_frags <- as.numeric(system(paste0(samtools, " view -f 0x40 -F 0x4 -F 0x8 ", test_file, " | cut -f9"), intern = TRUE))
prop_pairs <- as.numeric(system(paste0(samtools, " view -f 0x40 -f 0x2  ", test_file, " | cut -f9"), intern = TRUE))

all <- as.tibble(list(frags = abs(all_frags))) 
prop <- as.tibble(list(frags = abs(prop_pairs)))

all %>% 
  mutate(frags = abs(frags)) %>% 
  filter(frags < 1000 ) 

prop %>% 
  mutate(frags = abs(frags)) %>% 
  filter(frags < 1000) 

ggplot() +
  geom_density(data = all, aes(x = frags)) +
  geom_density(data = prop, aes(x = frags))
