#### annotate integrations imported by postprocess.R ####

## annotate the nearest feature in the specified bed files in array `nearest_bed`
## only include the distance to the nearest feature 

#### write integration bed file ####
# make necessary changes to ints tibble
bed <- ints %>% 
  dplyr::mutate(strand = ".") %>% 
  dplyr::select(Chr, IntStart, IntStop, strand, ReadID) %>%
  dplyr::arrange(Chr, IntStart)

# generate file names for saving
sorted_file <- paste0(tools::file_path_sans_ext(args[1]), ".filt.sorted.bed")

# save bed file
readr::write_tsv(bed, 
                 path = sorted_file,
                 col_names = FALSE)


#### loop over bed files and annotate nearest feature ####
for (bed in nearest_bed)
{
  
  cat("annotating integrations with nearest region from bed file", bed, "\n")
  
  # generate filename for sorted file
  nearest_bed <- paste0(tools::file_path_sans_ext(args[1]), ".nearest.", basename(tools::file_path_sans_ext(bed)), ".bed")

  
  # use bedtools to get nearest feature
  cmd <- paste0("bedtools closest -d -t first -a ", sorted_file, " -b ", bed, " > ", nearest_bed)
  system(cmd)
  
  # import nearest file to get nearest feature
  tmp <- readr::read_tsv(nearest_bed,
                         col_names = FALSE) 
  
  tmp_cols <- colnames(tmp)
  
  tmp <- tmp %>%
    dplyr::rename(Chr = X1) %>% 
    dplyr::rename(IntStart = X2) %>% 
    dplyr::rename(IntStop = X3) %>% 
    dplyr::rename(ReadID = X5) %>% 
    dplyr::rename_at(vars(starts_with(tmp_cols[length(tmp_cols)])), ~paste0(basename(tools::file_path_sans_ext(bed)), "_closest_dist")) %>% 
    dplyr::select(-starts_with("X"))  
  
  # join this annotation with ints
  ints <- left_join(ints, tmp, by = c("Chr", "IntStart", "IntStop", "ReadID"))
  
  # remove temp files
  file.remove(nearest_bed)
  
}

# remove temp bed files
file.remove(sorted_file)
