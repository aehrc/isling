#### Packages ####
library(tidyverse)

#import data
data <- read_tsv(file = "test_data/m366T.rearrange.txt") 

#make summary table of number of vector rearrangements
data %>% 
  group_by(possibleVecRearrange) %>% 
  count()

#make pie chart of possible vector rearrange
colns <- c("possibleVecRearrange", "Gaps")
for (i in colns) {
  count(data, !!ensym(i)) %>% 
    mutate( freq = n/sum(n) ) %>% 
    ggplot(aes(x="", y=freq, fill=!!ensym(i))) +
    geom_bar(stat = "identity", width=1)+
    geom_label(aes(label = paste0(n, "\n", round(freq*100), "%")), 
               position = position_fill(vjust = 0.5),
               show.legend = FALSE) +
    coord_polar("y", start=0) +
    theme_void() +
    theme(legend.position = "bottom")
  Sys.sleep(3)
} #end of loop over columns
  
  
#compare to FRG vec rearranges

soft <- read_tsv(file = "test_data/FRG206.soft.txt") 
count(soft, PossibleVectorRearrangement) %>% 
  mutate( freq = n/sum(n) ) %>% 
  ggplot(aes(x="", y=freq, fill=PossibleVectorRearrangement)) +
  geom_bar(stat = "identity", width=1)+
  geom_label(aes(label = paste0(n, "\n", round(freq*100), "%")), 
             position = position_fill(vjust = 0.5),
             show.legend = FALSE) +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.position = "bottom")


#make pie chart of number of gaps

  