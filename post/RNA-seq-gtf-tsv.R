#### annotate integrations imported by postprocess.R ####

## annotate the RNA-seq count of the nearest gene 
# using counts from a gtf file and counts from a specified column in a tsv file
## only include the "type" and "ID" columns from the gtf, 
## as well as the distance to this feature 

# make use of these vectors from the main script:
# gtf file to use: nearest_gtf_RNA
# name of RNA-seq column to use: RNA_seq_col
# RNA-seq tsv: RNA_seq_tsv

#### write integration bed file ####
# make necessary changes to ints tibble
bed <- ints %>% 
  dplyr::mutate(strand = ".") %>% 
  dplyr::select(Chr, IntStart, IntStop, strand, ReadID) %>%
  dplyr::arrange(Chr, IntStart)

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


#### loop over gtf, tsv files and annotate count of nearest feature ####

for (i in seq(length(nearest_gtf_RNA)))
{
  # get files and column for this 
  gtf <- nearest_gtf_RNA[i]
  tsv <- RNA_seq_tsv[i]
  col_name <- RNA_seq_col[i]
  
  cat("annotating RNA-seq counts from column", col_name, "in file", tsv, "with nearest region from gtf file", gtf, "\n")
  
  # check if we have already added the nearest gene
  if (paste0(basename(tools::file_path_sans_ext(gtf)), "_ID") %in% colnames(ints)) {
    
    # no need to add again
    tmp <- ints %>% 
      select(Chr, IntStart, IntStop, ReadID,
             paste0(basename(tools::file_path_sans_ext(gtf)), "_ID"), 
             paste0(basename(tools::file_path_sans_ext(gtf)), "_closest_dist"))
    
  }
  else {
    
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
  
    # remove temp files
    file.remove(nearest_gtf)
  }
  
  tmp <- tmp %>% 
    dplyr::mutate(gene_id = stringr::str_match(!!sym(paste0(basename(tools::file_path_sans_ext(gtf)), "_ID")), 
                                               "gene_id \\\"(\\w+)\\.\\d+\\\";")[,2]) %>% 
    dplyr::select(-contains(paste0(basename(tools::file_path_sans_ext(gtf)), "_ID")))
  
  # import RNA-seq tsv and add counts from specified column
  RNAseq_file <- tools::file_path_sans_ext(basename(tsv))
  RNAseq <- readr::read_tsv(tsv) %>% 
    setNames(paste0(RNAseq_file, "_", names(.))) 
  
  RNAseq <- RNAseq %>% 
    dplyr::mutate(gene_id = stringr::str_match(!!sym(paste0(RNAseq_file, "_gene_id")), '([\\w]+)\\.\\d+')[,2]) 
  
  RNAseq <- RNAseq %>% 
    dplyr::select(gene_id, !!sym(paste0(RNAseq_file, "_", col_name))) %>% 
    dplyr::filter(!is.na(gene_id))
  
  # annotate ints with RNA
  tmp <- tmp %>% 
    dplyr::left_join(RNAseq, by = "gene_id") %>% 
    dplyr::select(matches(paste("Chr", "IntStart", "IntStop", "ReadID", "gene_id", "dist", col_name, sep="|")))
  
  tmp <- tmp %>% 
    dplyr::rename(!!paste0(basename(tools::file_path_sans_ext(gtf)), "_gene_id") := gene_id) 
  
  # join this annotation with ints
  ints <- ints %>% 
    dplyr::left_join(tmp, by = c("Chr", "IntStart", "IntStop", "ReadID"))
  
  # remove duplicate column "<gtf file>_closest_dist
  ints <- ints %>% 
    dplyr::select(-!!sym(paste0(basename(tools::file_path_sans_ext(gtf)), "_closest_dist.y"))) 
  
  ints <- ints %>% 
    dplyr::rename(!!paste0(basename(tools::file_path_sans_ext(gtf)), "_closest_dist") := paste0(basename(tools::file_path_sans_ext(gtf)), "_closest_dist.x"))
  
  # remove duplicate column "<gtf file>_gene_id
  if (paste0(basename(tools::file_path_sans_ext(gtf)), "_gene_id.y") %in% colnames(ints)) {
    ints <- ints %>% 
      dplyr::select(-!!sym(paste0(basename(tools::file_path_sans_ext(gtf)), "_gene_id.y"))) 
  
    ints <- ints %>% 
      dplyr::rename(!!paste0(basename(tools::file_path_sans_ext(gtf)), "_gene_id") := paste0(basename(tools::file_path_sans_ext(gtf)), "_gene_id.x"))
  }
  # calculate count / dist
  count_col <- paste0(RNAseq_file, "_", col_name)
  dist_col <- paste0(basename(tools::file_path_sans_ext(gtf)), "_closest_dist")
  count_dist_col <- paste0(RNAseq_file, "_", col_name, "_divide_dist")
  ints <- ints %>% 
    dplyr::mutate(!!count_dist_col := ifelse(!!sym(dist_col) == 0, !!sym(count_col), !!sym(count_col)/!!sym(dist_col)))
  
}

# remove temp bed files
file.remove(sorted_file)
