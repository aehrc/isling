#### de-duplicate integrations imported by postprocess.R ####

# reatain only integrations with distinct read sequences
# consider only exact duplicates

# get counts of each read
counts <- ints %>% 
  dplyr::count(merged) %>% 
  dplyr::rename(n_duplicates = n)

# retain only distinct rows
ints <- ints %>% 
  dplyr::distinct(merged, .keep_all = TRUE)

# join counts info with ints
ints <- ints %>% 
  dplyr::left_join(counts, by="merged")

cat(nrow(ints), "integrations remaining after de-duplication based on read sequences\n")