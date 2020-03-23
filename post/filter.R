#### filter integrations imported by postprocess.R ####

# filtering criteria:
# not more than 20 ambiguous bases or are discordant read-pairs
# host edit distance not more than 5
# viral edit distance not more than 5
# not possible vector rearrangements
# not possilbe host translations

ints <- ints %>% 
  dplyr::filter(HostEditDist <= 5) %>% 
  dplyr::filter(ViralEditDist <= 5) %>% 
  dplyr::filter(NoAmbiguousBases <= 20 | OverlapType == 'discordant') %>% 
  dplyr::filter(PossibleVectorRearrangement == 'no') %>% 
  dplyr::filter(PossibleHostTranslocation == "no") 

cat(nrow(ints), "integrations remaining after filtering\n")