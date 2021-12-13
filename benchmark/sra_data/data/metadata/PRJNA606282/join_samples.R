library(tidyverse)

samples <- read_tsv("sampleDetails.tsv") %>%
 mutate(sampleName = str_replace(sampleName, "_.+$", "")) %>%
 mutate(sampleName = str_replace(sampleName, "^p", ""))

 
sra <- read_csv("SraRunTable.txt") %>%  
 mutate(`Library Name` = str_replace(`Library Name`, "_lib$", "")) %>%
 left_join(samples, by=c("Library Name"="sample"))

for (i in unique(sra$sampleName)) {
	print(glue::glue("animal id: {i}"))
	
	sra %>%
	 filter(sampleName == i) %>%
	 pull(Run) %>%
	 print()
}
