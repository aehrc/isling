#### annotate integrations imported by postprocess.R ####

## annotate the nearest feature in the specified gtf files in array `nearest_gtf`
## only include the "type" and "ID" columns from the gtf, 
## as well as the distance to this feature 

#### write integration bed file ####
# make necessary changes to ints tibble
bed <- ints %>% 
  dplyr::mutate(strand = ".") %>% 
  dplyr::select(Chr, IntStart, IntStop, strand, ReadID) 

# generate file names for saving
bed_file <- paste0(tools::file_path_sans_ext(args[1]), ".filt.bed")
sorted_file <- paste0(tools::file_path_sans_ext(args[1]), ".filt.sorted.bed")

# save bed file
readr::write_tsv(bed, 
                 path = bed_file,
                 col_names = FALSE)
                 
# sort bed file
cmd <- paste0("bedtools sort -i ", bed_file, " > ", sorted_file)
system(cmd)

#### loop over gtf files and annotate nearest feature ####
for (gtf in nearest_gtf)
{
  cat("annotating integrations with nearest region from gtf file", gtf, "\n")
  
  # generate filename for sorted file
  nearest_gtf <- paste0(tools::file_path_sans_ext(args[1]), ".nearest.", basename(tools::file_path_sans_ext(gtf)), ".bed")
  
  # use bedtools to get nearest feature
  cmd <- paste0("bedtools closest -d -t first -a ", sorted_file, " -b ", gtf, " > ", nearest_gtf)
  system(cmd)
  
  # import nearest file to get nearest feature
  bedCols <- c("Chr", "IntStart", "IntStop", "strand",  "ReadID")
  GTFcols <- c("discard1", "discard2", "discard3", "discard4", "discard5", "discard6", "discard7", "discard8", "ID", "closest_dist")
  colns <- c(bedCols, paste0(basename(tools::file_path_sans_ext(gtf)), "_", GTFcols))
  types = "ciiccccciicccci"
  tmp <- readr::read_tsv(nearest_gtf,
                  col_names = colns,
                  col_types = types) 
  tmp <- tmp %>% 
    dplyr::select(-contains("discard"))  %>% 
    dplyr::select(-strand)
  
  # join this annotation with ints
  ints <- left_join(ints, tmp, by = c("Chr", "IntStart", "IntStop", "ReadID"))
  
  # remove temp files
  file.remove(nearest_gtf)
  
}

# remove temp bed files
file.remove(sorted_file)
file.remove(bed_file)
