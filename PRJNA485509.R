

#### Packages ####
library(stringr)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(readr)
library(writexl)
library(circlize)
 

#### import data ####
raw_path = "../data/reads"
data_path = "../out"
out_path = "../out/summary/"

#get metadata from run tables
runTables <- list.files(raw_path, pattern ="SraRunTable.txt", recursive=TRUE)

#import metadata
metadata <- tibble(metaFile = runTables) %>%
  mutate(runTable = map(metaFile, ~ read_tsv(file.path(raw_path, .)))) %>%
  mutate(dataset = dirname(metaFile)) %>%
  unnest() %>%
  arrange(Run)

#get data from merged bed files
mergedbed_cols <- c('chr', 'start', 'stop', 'nReads', 'readNames')
merged <- metadata %>%
  mutate(beds = flatten_chr(map2(dataset, Run, ~ list.files(paste0(data_path, "/", .x, "/ints/", sep = ""), pattern=paste0(.y, ".+integrations.merged.bed", sep=""))))) %>%
  mutate(beds = flatten_chr(map2(dataset, beds, ~ paste0(data_path, "/", .x, "/ints/", .y, sep="")))) %>%
  mutate(merged = map(beds, ~read_tsv(., col_names = mergedbed_cols))) %>%
  mutate(uniq_sites = map_int(merged, ~ nrow(.)))

bed_cols <- c('chr', 'start', 'stop', 'readID')
bed <- metadata %>%
  mutate(beds = flatten_chr(map2(dataset, Run, ~ list.files(paste0(data_path, "/", .x, "/ints/", sep = ""), pattern=paste0(.y, ".+integrations.bed", sep=""))))) %>%
  mutate(beds = flatten_chr(map2(dataset, beds, ~ paste0(data_path, "/", .x, "/ints/", .y, sep="")))) %>%
  mutate(merged = map(beds, ~read_tsv(., col_names = bed_cols))) %>%
  mutate(uniq_sites = map_int(merged, ~ nrow(.)))

ints <- metadata %>%
  mutate(out = flatten_chr(map2(dataset, Run, ~ list.files(paste0(data_path, "/", .x, "/ints/", sep = ""), pattern=paste0(.y, ".+integrations.txt", sep=""))))) %>%
  mutate(out = flatten_chr(map2(dataset, out, ~ paste0(data_path, "/", .x, "/ints/", .y, sep="")))) %>%
  mutate(merged = map(out, ~read_tsv(.))) %>%
  mutate(uniq_sites = map_int(merged, ~ nrow(.)))

#save summary of data with number of sites included
for (i in unique(merged$dataset)) {
  toWrite <- list()
  data_filt <- merged %>% 
    filter(dataset == i)
  toWrite[[i]] <- data_filt  
}
write_xlsx(toWrite, path = paste(out_path, "summary_numSites.xlsx", sep = ""))

#plot number of sites by tissue
for (i in unique(merged$dataset)){
 merged %>%
  filter(dataset == i) %>%
  ggplot(aes(x=tissue, y = uniq_sites)) +
    geom_point(alpha = 0.5) +
    facet_wrap(vars(NucleicAcid)) +
    theme(axis.text.x = element_text(angle=45, hjust = 1)) 
  ggsave(paste(out_path, "summary_", i, "_mergedNumSites.pdf", sep = ""))
}

#plot only samples with less than 100 sites
for (i in unique(merged$dataset)) {
 merged %>%
  filter(dataset == i) %>%
  filter(uniq_sites < 100) %>%
    ggplot(aes(x = tissue, y = uniq_sites)) +
    geom_jitter(width = 0.1, height = 0, alpha = 0.5) + 
    facet_wrap(vars(NucleicAcid)) +
    theme(axis.text.x = element_text(angle=45, hjust = 1)) 
  ggsave(paste(out_path, "summary_", i, "_mergedNumSites100orLess.pdf", sep = ""))
}


###unnest tibbles to plot individual sites ####

merged <- unnest(merged)
bed <- unnest(bed)
ints <- unnest(ints)

#plot location of merged bed sites for each tissue
for (i in unique(merged$dataset)) {
 filt <- filter(merged, dataset == i)
  host <- filt %>% 
    select(Organism) %>% 
    unique()
  if (str_detect(host, "Mus musculus")) {host <- "mm10"}
  chroms <- if (str_detect(host, "hg38")) {paste0("chr", c(1:22, "X", "Y"))} else {paste0("chr", c(1:20, "X", "Y"))}
  
  for (k in unique(filt$NucleicAcid)) {
    filt2 <- filter(filt, NucleicAcid == k)
    for (j in unique(filt2$tissue)) {
      bed <- filter(filt2, tissue == j) %>%
        select(chr, start, stop) %>%
        as.data.frame()
        
      pdf(paste0(out_path, "hostCircos_", i, "_", k, "_", j, ".pdf", sep = ""))
      if (str_detect(host, "hg38")) {circos.initializeWithIdeogram(species = host, 
                                                                   chromosome.index = chroms)} 
      else {circos.initializeWithIdeogram(species = host)}
      circos.genomicRainfall(bed, track.height = 0.2)
      circos.genomicDensity(bed, track.height = 0.1)
      dev.off()
      circos.clear()
    }
  }
}




