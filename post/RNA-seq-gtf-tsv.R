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
sorted_file <- paste0(tools::file_path_sans_ext(args[1]), ".filt.sorted.bed")

# save bed file
readr::write_tsv(bed, 
                 path = sorted_file,
                 col_names = FALSE)

#### loop over gtf, tsv files and annotate count of nearest feature ####

for (gtf in nearest_gtf)
{
  gtf <- nearest_gtf_RNA[i]
  tsv <- RNA_seq_tsv[i]
  col_name <- RNA_seq_col[i]
  
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
    dplyr::select(-strand) %>% 
    dplyr::mutate(gene_id = stringr::str_match(ID, 'gene_id \\"([\\w]+)\\.\\d+\\";')[,2]) %>% 
    dplyr::select(-ID)
  
  # import RNA-seq tsv and add counts from specified column
  RNAseq <- readr::read_tsv(tsv) %>% 
    dplyr::mutate(gene_id = str_match(gene_id, '([\\w]+)\\.\\d+')[,2]) %>% 
    dplyr::filter(!is.na(gene_id))
  
  #rename columns to add filename
  RNAseq_file <- basename(tsv)
  RNAseq <- RNAseq %>% 
    setNames(paste0(RNAseq_file, "_", names(.)))
  
  # annotate ints with RNA
  tmp <- tmp %>% 
    dplyr::left_join(tmp, RNAseq, by = "gene_id") %>% 
    dplyr::select(matches(paste("Chr", "IntStart", "IntStop", "ReadID", "gene_id", "dist", col_name, sep="|")))
  
  # join this annotation with ints
  ints <- dplyr::left_join(ints, tmp, by = c("Chr", "IntStart", "IntStop", "ReadID"))
  
  # remove temp files
  file.remove(nearest_gtf)
  
}

# remove temp bed files
file.remove(sorted_file)