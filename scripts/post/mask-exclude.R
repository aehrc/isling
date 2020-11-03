#### filter integrations imported by postprocess.R ####

## mask according to bed file(s) in array `mask-exclude`
## exclude any integrations that overlap with the regions in this these files

# function to check for overlaps
overlaps <- function(chr, start, stop, rchr, rstart, rstop)
{
  overlap <- dplyr::case_when(
  chr != rchr ~ FALSE,
  (stop >= rstart) & (stop <= rstop) ~ TRUE, # stop is between rstart and rstop
  (start >= rstart) & (start <= rstop) ~ TRUE, # start is between rstart and rstop
  (start <= rstart) & (stop >= rstop) ~ TRUE, # start is after rstart and stop is after rstart
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
  if (sum(unique(bed_regions$X1) %in% unique(ints$Chr)) == 0)
  {
    cat("no overlap between chromosome names in ints file and ", bed)
  }
  
  # for each line, filter integrations that are in this region
  for (i in 1:nrow(bed_regions))
  {
    chr <- bed_regions$X1[i]
    start <- bed_regions$X1[i]
    stop <- bed_regions$X1[i]
    
    ints <- ints %>% 
      dplyr::filter(!overlaps(Chr, IntStart, IntStop, chr, start, stop))
    
  }
  cat(nrow(ints), "integrations left after excluding integrations in regions specified in file", bed, "...\n")
}



#testchr <- c("chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr2")
#teststart <- c(1, 4, 3, 4, 1, 0, 7, 5)
#teststop <- c(3, 5, 4, 6, 6, 1, 8, 6)
#overlaps(testchr, teststart, teststop, "chr1", 2, 5)
