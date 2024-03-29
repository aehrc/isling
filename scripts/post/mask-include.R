#### filter integrations imported by postprocess.R ####

## mask according to bed file(s) in array `mask-exclude`
## only include integrations that overlap with the regions in this these files

# function to check for overlaps
overlaps <- function(chr, start, stop, rchr, rstart, rstop)
{
  overlap <- dplyr::case_when(
    chr != rchr ~ FALSE,
    (stop >= rstart) & (stop <= rstop) ~ TRUE,
    (start >= rstart) & (start <= rstop) ~ TRUE,
    (start <= rstart) & (stop >= rstop) ~ TRUE,
    TRUE ~ FALSE
  )
  return(overlap)
}

# loop through bed files and mask each one
for (bed in mask_exclude)
{
  # import regions from bed file
  bed_regions <- readr::read_tsv(file = bed,
                          col_names = FALSE,
                          col_types = cols(
                            X1 = col_character(),
                            X2 = col_double(),
                            X3 = col_double())
                          )
  
  # check that there is some overlap between chromosome names in bed_regions and ints
  if (sum(bed_regions$X1 %in% ints$Chr))
  {
    stop(paste0("no overlap between chromosome names in ints file and ", bed))
  }
  
  # for each line, filter integrations that are in this region
  for (i in 1:nrow(bed_regions))
  {
    chr <- bed_regions$X1[i]
    start <- bed_regions$X1[i]
    stop <- bed_regions$X1[i]
    
    ints <- ints %>% 
      dplyr::filter(overlaps(Chr, IntStart, IntStop, chr, start, stop))
    
  }
  cat(nrow(ints), "integrations left after excluding integrations in regions specified in file", bed, "...\n")
}



#testchr <- c("chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr2")
#teststart <- c(1, 4, 3, 4, 1, 0, 7, 5)
#teststop <- c(3, 5, 4, 6, 6, 1, 8, 6)
#overlaps(testchr, teststart, teststop, "chr1", 2, 5)
