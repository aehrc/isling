#### de-duplicate integrations imported by postprocess.R ####

# reatain only integrations with distinct read sequences
# consider only exact duplicates

ints <- ints %>% 
  dplyr::distinct(merged, .keep_all = TRUE)

cat(nrow(ints), "integrations remaining after de-duplication based on read sequences\n")