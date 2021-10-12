## helper functions for summary Rmd

max_sample_length <- 10
max_contigs_plot <- 25

pal <- wesanderson::wes_palette("Royal1")

col_spec <- readr::cols(
  .default = readr::col_character(),
  IntStart = readr::col_double(),
  IntStop = readr::col_double(),
  VirusStart = readr::col_double(),
  VirusStop = readr::col_double(),
  NoAmbiguousBases = readr::col_double(),
  HostEditDist = readr::col_double(),
  ViralEditDist = readr::col_double(),
  TotalEditDist = readr::col_double(),
  PossibleHostTranslocation = readr::col_logical(),
  PossibleVectorRearrangement = readr::col_logical(),
  HostAmbiguousLocation = readr::col_logical(),
  ViralAmbiguousLocation = readr::col_logical(),
  HostMapQ = readr::col_double(),
  ViralMapQ = readr::col_double()
)


# https://stackoverflow.com/questions/61711777/how-to-create-an-adaptative-list-ulli-using-htmltools
vec_to_html_list <- function(vec) {
  return(
    htmltools::tags$ul(
      purrr::map(vec, ~htmltools::tags$li(htmltools::tags$code(.)))
    )
  )
}

get_contig_lengths <- function(bwa_index_prefix) {
  filepath <- paste(bwa_index_prefix, ".ann", sep="")
  con = file(filepath, "r")
  i <- 1
  contigs <- list()
  # first line contains summary for whole genome
  readLines(con, n = 1)
  while ( TRUE ) {
    
    # get next line
    line <- readLines(con, n = 1)
    
    # break if we've reached the end of the file
    if ( length(line) == 0 ) {
      break
    }   
    
    # split line into fields
    line <- stringr::str_split(line, "\\s", simplify=TRUE)[1,]
    
    if (i %% 2 == 1) {
      # get name of contig
      contig <- line[2]
    } else {
      # get length of contig from next line
      contigs[[contig]] = as.double(line[2])
    }
    
    i <- i + 1
    
  }
  
  close(con)
  return(contigs)
  
}



get_top_contigs <- function(ints, col_name, contig_lengths) {
  contigs_join <- tibble::tibble(
    name = names(contig_lengths),
    length = unlist(unname(contig_lengths))
  )
  
  right_join <- "name"
  
  include_group <- ifelse(col_name == "Chr", "ambig. (virus)", "ambig. (host)")
  
  contigs_plot <- ints %>% 
    dplyr::filter(type == "unique" | type == include_group) 
  
  if (nrow(contigs_plot) == 0) {
    return(contigs_plot)
  }
  
 contigs_plot <- contigs_plot %>% 
    dplyr::group_by(!!dplyr::sym(col_name)) %>% 
    dplyr::summarise(n = dplyr::n()) %>% 
    dplyr::arrange(desc(n)) %>% 
    head(max_contigs_plot) 
 
 if (nrow(contigs_plot) == 0) {
   return(contigs_plot)
 }
 
 contigs_plot <- contigs_plot %>% 
    dplyr::left_join(contigs_join, by=setNames(right_join, col_name)) %>% 
    dplyr::rename(contig_name = !!dplyr::sym(col_name)) %>% 
    dplyr::mutate(start=1) %>% 
    dplyr::arrange(factor(contig_name, levels=stringr::str_sort(contig_name, numeric=TRUE)))
  
  return(contigs_plot)
}
