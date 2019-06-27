#### Plot integration site FRG206 ####

#### Packages ####
library(tidyverse)
library(circlize)
  
#working with just one dataset  
FRG206 <- read_tsv(file = "data/FRG206.txt", col_names = TRUE)

#add a few extra summary columns about ambiguity
FRG206_extra <- FRG206 %>% 
  mutate(has_issue = ifelse((PossibleVectorRearrangement == "yes" | 
                        PossibleHostTranslocation == "yes" |
                          HostPossibleAmbiguous == "yes" |
                          ViralPossibleAmbiguous == "yes"), 'yes', 'no')) %>% 
  mutate(has_rearrangement = ifelse(PossibleVectorRearrangement == "yes" | 
                                      PossibleHostTranslocation == "yes", "yes", "no")) %>% 
  mutate(has_locAmbig = ifelse( HostPossibleAmbiguous == "yes" |
                                  ViralPossibleAmbiguous == "yes", 'yes', 'no')) 

#plot frequency of yes/no answers for each possible ambiguity
colns <- c("PossibleVectorRearrangement", "PossibleHostTranslocation", "HostPossibleAmbiguous", "ViralPossibleAmbiguous", "has_issue", "has_rearrangement", "has_locAmbig")
for (i in colns) {
  count(FRG206_extra, !!ensym(i)) %>% 
    mutate( freq = n/sum(n) ) %>% 
    ggplot(aes(x="", y=freq, fill=!!ensym(i))) +
    geom_bar(stat = "identity", width=1)+
    geom_label(aes(label = paste0(n, "\n", round(freq*100), "%")), 
               position = position_fill(vjust = 0.5),
               show.legend = FALSE) +
    coord_polar("y", start=0) +
    theme_void() +
    theme(legend.position = "bottom")
  ggsave(paste0("out/",gsub(" ", "", i), ".pdf"))
} #end of loop over columns

#### check out level of ambiguity in FRG206 ####
#rearrange FRG206 for plotting
count_ambig <- FRG206 %>% 
  select(PossibleVectorRearrangement, PossibleHostTranslocation, HostPossibleAmbiguous, ViralPossibleAmbiguous) %>% 
  gather(key = "type", value = "present") %>% 
  mutate(host_or_vector = ifelse(str_detect(type, "Host"), "host", "virus")) %>% 
  mutate(type = recode(type, 
                       "PossibleVectorRearrangement" = "rearrange",
                       "PossibleHostTranslocation" = "rearrange",
                       "HostPossibleAmbiguous" = "ambig_loc",
                       "ViralPossibleAmbiguous" = "ambig_loc")) 
  
#make plot
count_ambig %>% 
    group_by(type, host_or_vector, present) %>% 
    summarise( counts = n()) %>%
    ggplot(aes(x = "", y = counts, fill = present)) +
    geom_bar(stat = "identity") +
    coord_polar("y", start = 0) +
    facet_grid(rows = vars(host_or_vector), cols = vars(type)) 
ggsave("out/ambig_facet.pdf")


#### how is ambiguity distributed across human genome? ####
FRG206_extra %>% 
  group_by(Chr, has_issue) %>% 
  summarise(counts = n()) %>% 
  ggplot(aes (x = Chr, y = counts, fill = has_issue)) +
  geom_bar(stat = "identity", position = "dodge") + 
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=90, hjust=1)) 

#get bed of flagged sites
bed_issues <- FRG206_extra %>% 
  filter(has_issue == "yes") %>% 
  select(1:3) %>% 
  filter(str_detect(Chr, paste0("chr", c(1:22, "X", "Y"), "$", collapse = "|"))) %>% 
  as.FRG206.frame()

#get bed of sites with no flag
bed_noissues <- FRG206_extra %>% 
  filter(has_issue == "no") %>% 
  select(1:3) %>% 
  filter(str_detect(Chr, paste0("chr", c(1:22, "X", "Y"), "$", collapse = "|"))) %>% 
  as.FRG206.frame()

bed_list <- list(bed_issues, bed_noissues)

#make plot
pdf("out/density.pdf")
circos.initializeWithIdeogram(species='hg38', chromosome.index = paste0("chr", c(1:22, "X", "Y")))
circos.genomicRainfall(bed_list, col = c("#00BFC4", "#F8766D"), track.height = 0.2)
circos.genomicDensity(bed_issues, col = c("#00BFC4"), track.height = 0.1)
circos.genomicDensity(bed_noissues, col = c("#F8766D"), track.height = 0.1)
dev.off()
circos.clear()

#exclude chromosome 14 so can get a sense of rest of the sites
bed_issues_no14 <- FRG206_extra %>% 
  filter(!str_detect(Chr, "chr14")) %>% 
  filter(has_issue == "yes") %>% 
  select(1:3) %>% 
  filter(str_detect(Chr, paste0("chr", c(1:22, "X", "Y"), "$", collapse = "|"))) %>% 
  as.FRG206.frame()

#get bed of sites with no flag
bed_noissues_no14 <- FRG206_extra %>% 
  filter(!str_detect(Chr, "chr14")) %>% 
  filter(has_issue == "no") %>% 
  select(1:3) %>% 
  filter(str_detect(Chr, paste0("chr", c(1:22, "X", "Y"), "$", collapse = "|"))) %>% 
  as.FRG206.frame()

bed_list2 = list(bed_issues_no14, bed_noissues_no14)
pdf("out/density_no14.pdf")
circos.initializeWithIdeogram(species='hg38', chromosome.index = paste0("chr", c(1:22, "X", "Y")))
circos.genomicRainfall(bed_list2, col = c("#00BFC4", "#F8766D"), track.height = 0.2)
circos.genomicDensity(bed_issues_no14, col = c("#00BFC4"), track.height = 0.1)
circos.genomicDensity(bed_noissues_no14, col = c("#F8766D"), track.height = 0.1)
dev.off()
circos.clear()


#make plot with links between virus and human chromosomes
df = FRG206.frame(
  name  = c("TP53",  "TP63",    "TP73"),
  start = c(7565097, 189349205, 3569084),
  end   = c(7590856, 189615068, 3652765))
circos.genomicInitialize(df)
circos.initializeWithIdeogram(species='hg38', chromosome.index = paste0("chr", c(1:22, "X", "Y")))



#get non-ambiguous FRG206
not_ambig <- FRG206 %>% 
  filter(PossibleVectorRearrangement == "no") %>% 
  filter(PossibleHostTranslocation == "no") %>% 
  filter(HostPossibleAmbiguous == "no") %>% 
  filter(ViralPossibleAmbiguous == "no") 




####working with all datasets ####
data_path = "data/2019-06-26/"
out_path = "out/2019-06-26/"
files <- dir(data_path, pattern ="integrations.txt")

#import all datasets
data <- data_frame(filename = files) %>% # create a data frame holding the file names
  mutate(integrations = map(filename, ~ read_tsv(file.path(data_path, .), na = c("", "NA", "?"), col_types = cols(Chr = col_character(), .default = col_guess()))) # a new data column
  ) 

data

#add column with number of integration sites for each dataset
data <- data %>% 
  mutate( n = sapply(data$integrations, nrow)) %>% 
  mutate( dataset = ifelse(str_detect(filename, "FRG"), "FRG", 
                           ifelse(str_detect(filename, "combinedReads"), "patient", 
                                  ifelse(str_detect(filename, "m3"), "mouse", "second_patient"))))

#add column with host used for alignment for later use
data <- data %>% 
  mutate( host = ifelse(str_detect(dataset, "mouse"), "mm10", "hg38"))


#which datasets have no sites?
no_sites <- filter(data, n == 0)
no_sites

#retain only datasets with some integration sites
data <- data %>% 
  filter(n > 0)  %>% 
  unnest() %>% 
  mutate(has_issue = ifelse((PossibleVectorRearrangement == "yes" | 
                               PossibleHostTranslocation == "yes" |
                               HostPossibleAmbiguous == "yes" |
                               ViralPossibleAmbiguous == "yes"), 'yes', 'no')) %>% 
  mutate(has_rearrangement = ifelse(PossibleVectorRearrangement == "yes" | 
                                      PossibleHostTranslocation == "yes", "yes", "no")) %>% 
  mutate(has_locAmbig = ifelse( HostPossibleAmbiguous == "yes" |
                                  ViralPossibleAmbiguous == "yes", 'yes', 'no'))

#get 'unique' integrations
data <- data %>% distinct(filename, Chr, IntStart, IntStop, VirusRef, VirusStart, VirusStop, NoAmbiguousBases, OverlapType, Orientation, PossibleHostTranslocation, PossibleVectorRearrangement, HostPossibleAmbiguous, ViralPossibleAmbiguous,
                          .keep_all=TRUE)

#get summary of number of sites per 
data %>% 
  group_by(dataset, filename) %>% 
  count()  %>% 
  write_tsv(paste0(out_path, "_summary/num_sites.txt"), na = "NA",
            quote_escape = "double")

#loop over columns to make pie charts of all columns
colns <- c("OverlapType", "PossibleVectorRearrangement", "PossibleHostTranslocation", "HostPossibleAmbiguous", "ViralPossibleAmbiguous", "has_issue", "has_rearrangement", "has_locAmbig")
for (i in colns) {
  data %>%
    filter(!str_detect(dataset, "second_patient")) %>% #ignore this dataset for now
    mutate(dataset = as.factor(dataset)) %>% #need to make this a factor for complete to work
    count(!!ensym(i), dataset) %>% # get counts by coln[i] and dataset
    complete(!!ensym(i), dataset, fill = list(n = 0)) %>% #fill in zero rows
    group_by(dataset) %>% #need this to get frequency of counts by dataset
      mutate( freq = n/sum(n) ) %>% #get frequency
      ggplot(aes(x="", y=freq, fill=!!ensym(i))) +
      geom_bar(stat = "identity", width=1) +
      geom_label(aes(label = paste0(n, "\n", round(freq*100), "%")), 
                 position = position_fill(vjust = 0.5),
                 show.legend = FALSE) +
      coord_polar("y", start=0) +
      theme_void() +
      theme(legend.position = "bottom") +
      facet_wrap(vars(dataset))
  ggsave(paste0("out/pie_",gsub(" ", "", i), ".pdf"))
}


#make circos plots for all datasets
datasets <- unique(data$dataset)
for (i in datasets) {
  
  #get host for this dataset
  host <- data %>% 
    filter(str_detect(dataset, i)) %>% 
    select(host) %>% 
    unique()
  
  #get chromosomes for host
  chroms <- if (str_detect(host, "hg38")) {paste0("chr", c(1:22, "X", "Y"))} else {paste0("chr", c(1:20, "X", "Y"))}
  
  #make bed file of data with no issues
  bed_noissues <- data %>% 
    filter(str_detect(dataset, i)) %>% 
    filter(has_issue == "no") %>% 
    select(Chr:IntStop) %>%
    mutate(Chr = paste0("chr", Chr)) %>% 
    filter(str_detect(Chr, paste(chroms, collapse = "|"))) %>% 
    as.data.frame()
  
  #make bed file of data with issues
  bed_issues <- data %>% 
    filter(str_detect(dataset, i)) %>% 
    filter(has_issue == "yes") %>% 
    select(Chr:IntStop) %>%
    mutate(Chr = paste0("chr", Chr)) %>% 
    filter(str_detect(Chr, paste(chroms, collapse = "|"))) %>% 
    as.data.frame()
  
  #make bed file of all data
  bed <- data %>% 
    filter(str_detect(dataset, i)) %>% 
    select(Chr:IntStop) %>%
    mutate(Chr = paste0("chr", Chr)) %>% 
    filter(str_detect(Chr, paste(chroms, collapse = "|"))) %>% 
    as.data.frame()
  
  #make plot of all data
  pdf(paste0("out/hostCircos_", i, "_density_all.pdf"))
  if (str_detect(host, "hg38")) {circos.initializeWithIdeogram(species = host, 
                                                               chromosome.index = chroms)} else {circos.initializeWithIdeogram(species = host)}
  circos.genomicRainfall(bed, track.height = 0.2)
  circos.genomicDensity(bed, track.height = 0.1)
  dev.off()
  circos.clear()
  
  #make plot of data with no issues
  bed_list2 = list(bed_issues, bed_noissues)
  pdf(paste0("out/hostCircos_", i, "_density_issues.pdf"))
  if (str_detect(host, "hg38")) {circos.initializeWithIdeogram(species = host, 
                                                               chromosome.index = chroms)}
  else {circos.initializeWithIdeogram(species = host)}
  circos.genomicRainfall(bed_list2, col = c("#00BFC4", "#F8766D"), track.height = 0.2)
  circos.genomicDensity(bed_issues, col = c("#00BFC4"), track.height = 0.1)
  circos.genomicDensity(bed_noissues, col = c("#F8766D"), track.height = 0.1)
  dev.off()
  circos.clear()

}









