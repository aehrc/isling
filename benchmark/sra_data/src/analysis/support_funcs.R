# scripts for common tasks for analysing external data

# import isling results
importIslingMergedSRA <- function(isling_dir, meta) {
  
  isling_suffix <- "merged.txt"
  
  coltypes <- cols(
    Chr = col_character(),
    IntStart = col_double(),
    IntStop = col_double(),
    Virus = col_character(),
    VirusStart = col_double(),
    VirusStop = col_double(),
    nChimeric = col_double(),
    nDiscordant = col_double(),
    SiteID = col_double(),
    ReadIDs = col_character()
  )
  
  isling_df <- tibble(
    f = list.files(isling_dir, isling_suffix, recursive=TRUE),
    data = map(f, ~read_tsv(file.path(isling_dir,.), col_types = coltypes))
  )  %>% 
    mutate(Run = str_split(basename(f), "\\.", simplify=TRUE)[,1]) %>% 
    mutate(host = str_split(basename(f), "\\.", simplify=TRUE)[,2]) %>%   
    mutate(virus = str_split(basename(f), "\\.", simplify=TRUE)[,3])  %>% 
    mutate(params = basename(dirname(dirname(f)))) %>% 
    left_join( meta, by="Run") %>% 
    unnest(data)
  
  return(isling_df)
  
}

importIslingVirusAmbig <- function(isling_dir, meta) {
  
  isling_suffix <- "virus_ambig.txt"
  
  coltypes <- cols(
    .default = col_character(),
    IntStart = col_double(),
    IntStop = col_double(),
    VirusStart = col_double(),
    VirusStop = col_double(),
    NoAmbiguousBases = col_double(),
    HostEditDist = col_double(),
    ViralEditDist = col_double(),
    TotalEditDist = col_double(),
    PossibleHostTranslocation = col_logical(),
    PossibleVectorRearrangement = col_logical(),
    HostAmbiguousLocation = col_logical(),
    ViralAmbiguousLocation = col_logical(),
    HostMapQ = col_double(),
    ViralMapQ = col_double()
  )
  
  isling_df <- tibble(
    f = list.files(isling_dir, isling_suffix, recursive=TRUE),
    data = map(f, ~read_tsv(file.path(isling_dir,.), col_types = coltypes, na = c("", "NA", "None"))))  %>% 
    mutate(Run = str_split(basename(f), "\\.", simplify=TRUE)[,1]) %>% 
    mutate(host = str_split(basename(f), "\\.", simplify=TRUE)[,2]) %>%   
    mutate(virus = str_split(basename(f), "\\.", simplify=TRUE)[,3])  %>% 
    mutate(params = basename(dirname(dirname(f)))) %>% 
    left_join( meta, by="Run") %>% 
    unnest(data)
  
  
  return(isling_df)
  
}

importIslingHostAmbig <- function(isling_dir, meta) {
  
  isling_suffix <- "host_ambig.txt|both_ambig.txt"
  
  coltypes <- cols(
    .default = col_character(),
    IntStart = col_double(),
    IntStop = col_double(),
    VirusStart = col_double(),
    VirusStop = col_double(),
    NoAmbiguousBases = col_double(),
    HostEditDist = col_double(),
    ViralEditDist = col_double(),
    TotalEditDist = col_double(),
    PossibleHostTranslocation = col_logical(),
    PossibleVectorRearrangement = col_logical(),
    HostAmbiguousLocation = col_logical(),
    ViralAmbiguousLocation = col_logical(),
    HostMapQ = col_double(),
    ViralMapQ = col_double()
  )
  
  isling_df <- tibble(
    f = list.files(isling_dir, isling_suffix, recursive=TRUE),
    data = map(f, ~read_tsv(file.path(isling_dir,.), col_types = coltypes, na = c("", "NA", "None"))))  %>% 
    mutate(Run = str_split(basename(f), "\\.", simplify=TRUE)[,1]) %>% 
    mutate(host = str_split(basename(f), "\\.", simplify=TRUE)[,2]) %>%   
    mutate(virus = str_split(basename(f), "\\.", simplify=TRUE)[,3])  %>% 
    mutate(params = basename(dirname(dirname(f)))) %>% 
    left_join( meta, by="Run") %>% 
    unnest(data)
  
  
  
  return(isling_df)
  
}

importIslingHostSRA <- function(isling_dir, meta) {
  
  merged_unique <- importIslingMergedSRA(isling_dir, meta)
  
  host_unique <- importIslingVirusAmbig(isling_dir, meta)
  
  host_unique <- host_unique %>% 
    select(Chr, IntStart, IntStop, Orientation, VirusRef, VirusStart, VirusStop,
           VirusOrientation, Type, ReadID, Run, host, virus, params) %>% 
    rename(Virus = VirusRef) %>% 
    mutate(nChimeric = ifelse(Type == "chimeric", 1, 0)) %>% 
    mutate(nDiscordant = ifelse(Type == "discordant", 1, 0)) %>% 
    rename(ReadIDs = ReadID) %>% 
    rename(JunctionType = Orientation) %>% 
    rename(VirusOri = VirusOrientation) %>% 
    select(-Type) %>% 
    mutate(SiteID = as.double(NA)) %>% 
    left_join(meta, by="Run")
    
  return(bind_rows(merged_unique, host_unique))
  
}

importIslingAllSRA <- function(isling_dir, meta) {
  # import unique integrations
  unique_host <- importIslingHostSRA(isling_dir, meta)
  # import ambiguous integrations
  ambig_host <- importIslingHostAmbig(isling_dir, meta)
  
  # unpack 'AltInts' column into tibble with all locations
  ambig_host_all <- tibble()
  
  for (i in seq(nrow(ambig_host))) {
    # get info from 'primary' integration site
    params <- c(ambig_host$params[[i]])
    runs <- c(ambig_host$Run[[i]])
    chrs <- c(ambig_host$Chr[[i]])
    starts <- c(ambig_host$IntStart[[i]])
    stops <- c(ambig_host$IntStop[[i]])
    oris <- c(ambig_host$Orientation[[i]])
    viruses <- c(ambig_host$VirusRef[[i]])
    virus_starts <- c(ambig_host$VirusStart[[i]])
    virus_stops <- c(ambig_host$VirusStop[[i]])
    virus_oris <- c( ambig_host$VirusOrientation[[i]])
    read_ids <- c(ambig_host$ReadID[[i]])
    
    # get alternate locations
    if (!(is.na(ambig_host$AltLocs[[i]]))) {
      for (int in str_split(ambig_host$AltLocs[[i]], ";")) {
        int <- str_split(int, ",")[[1]]
        
        # get info about location in host
        host <- str_split(int[1], ":")[[1]]
        chrs <- c(chrs, host[1])
        host <- str_split(host[2], "-")[[1]]
        starts <- c(starts, as.numeric(host[1]))
        stops <- c(stops, as.numeric(host[2]))  
        oris <- c(oris, int[2])
        
        # get info about location in virus
        virus <- str_split(int[3], ":")[[1]]
        viruses <- c(viruses, virus[1])
        virus <- str_split(virus[2], "-")[[1]]
        virus_starts <- c(virus_starts, as.numeric(virus[1]))
        virus_stops <- c(virus_stops , as.numeric(virus[2]))
        virus_oris <- c(virus_oris, int[4])
        
        #other stuff needed
        read_ids <- c(read_ids, ambig_host$ReadID[[i]])
        runs <- c(runs, ambig_host$Run[[i]])
        params <- c(params, ambig_host$params[[i]])
        
      }
    }
    
    # add current row to all dataframe
    ambig_host_all <- bind_rows(ambig_host_all, tibble(
      Chr = chrs,
      IntStart = starts,
      IntStop = stops,
      ori = oris,
      virus = viruses,
      virus_start = virus_starts,
      virus_stop = virus_stops,
      virus_ori = virus_oris,
      ReadIDs = read_ids,
      Run = runs,
      params = params
    )
    )
    
  }  
  
  # combine ambiugous and unique
  isling_all <- unique_host %>% 
    select(Chr, IntStart, IntStop, Run, ReadIDs, params) %>% 
    bind_rows(ambig_host) %>% 
    ungroup() %>% 
    mutate(params = paste0(params, "_unique-plus-ambig")) %>% 
    group_by(Run, ReadIDs, params) %>% 
    distinct(Chr, IntStart, IntStop) %>% 
    ungroup()
  
  return(left_join(isling_all, meta, by='Run'))
}




importPolyidus <- function(polyidus_dir) {
  polyidus_filename <- "exactHpvIntegrations.bed"
  
  coltypes <- cols(
    chrom = col_character(),
    chromStart = col_double(),
    chromEnd = col_double(),
    name = col_character(),
    score = col_double(),
    strand = col_character(),
    virusStart = col_double(),
    virusEnd = col_double(),
    virusStrand = col_character(),
    numReads = col_double(),
    bioSample = col_character()
  )
  
  polyidus <- tibble(
    f = list.files(polyidus_dir, pattern = polyidus_filename, recursive = TRUE),
    data = map(f, ~read_tsv(file.path(polyidus_dir, .), col_types = coltypes))
  ) %>% 
    mutate(info = dirname(dirname(f))) %>% 
    mutate(host = str_split(info, "\\.", simplify=TRUE)[,1]) %>% 
    mutate(virus = str_split(info, "\\.", simplify=TRUE)[,2]) %>% 
    mutate(Run = str_split(info, "\\.", simplify=TRUE)[,3])  %>% 
    mutate(params = "polyidus") %>% 
    unnest(data)
  
  return(polyidus)
  
}

# import polyidus results 
importPolyidusSRA <- function(polyidus_dir, meta) {
  
  polyidus <- importPolyidus(polyidus_dir) %>% 
    left_join(meta, by="Run") 
  
  return(polyidus)
  
}

importVifi <- function(vifi_dir) {
  vifi_filename <-"output.clusters.txt.range"
  
  column_types <- cols(
    Chr = col_character(),
    Min = col_double(),
    Max = col_double(),
    Split1 = col_double(),
    Split2 = col_double()
  )
  
  vifi <- tibble(
    f = file.path(vifi_dir, list.files(vifi_dir, pattern = vifi_filename, recursive = TRUE)),
    data = map(f, ~read_csv(., col_types = column_types)),
    Run = str_split(basename(dirname(f)), "\\.", simplify=TRUE)[,1],
    params = "vifi"
  ) %>% 
    unnest(data)
  
  return(vifi) 
}

#import vifi results
importVifiSRA <- function(vifi_dir, meta) {
  
  
  vifi <- importVifi(vifi_dir)
  
  return(left_join(vifi, meta, by="Run"))
}

importSeeksv <- function(seeksv_dir, host_fai, remove_chr=FALSE) {
  
  seeksv_pattern <- ".integrations.txt"
  
  files <- list.files(seeksv_dir, pattern = seeksv_pattern)
  
  print(glue::glue("importing {length(files)} files"))
  
  coltypes <- readr::cols(
    `@left_chr` = col_character(),
    left_pos = col_double(),
    left_strand = col_character(),
    left_clip_read_NO = col_double(),
    right_chr = col_character(),
    right_pos = col_double(),
    right_strand = col_character(),
    right_clip_read_NO = col_double(),
    microhomology_length = col_double(),
    abnormal_readpair_NO = col_double(),
    svtype = col_character(),
    left_pos_depth = col_double(),
    right_pos_depth = col_double(),
    average_depth_of_left_pos_5end = col_double(),
    average_depth_of_left_pos_3end = col_double(),
    average_depth_of_right_pos_5end = col_double(),
    average_depth_of_right_pos_3end = col_double(),
    left_pos_clip_percentage = col_double(),
    right_pos_clip_percentage = col_double(),
    left_seq_cigar = col_character(),
    right_seq_cigar = col_character(),
    left_seq = col_character(),
    right_seq = col_character()
  )
  
  # need to filter out anything that doesn't involve a host chromsome and a virus
  colnames <- c("chr", "length", "byte1", "byte2", "byte3")
  host_chrs <- readr::read_tsv(host_fai, col_names = colnames) %>% 
    dplyr::pull(chr) 
  host_chrs <- paste0("^", host_chrs, "$", collapse = "|")
  
  # import seeksv data and filter for integrations
  seeksv <- tibble::tibble(
    f = file.path(seeksv_dir, files),
    Run = str_split(basename(f), "\\.", simplify=TRUE)[,1],
    params = "seeksv",
    data = map(f, ~read_tsv(., col_types = coltypes))
  ) %>% 
    tidyr::unnest(data)  

  if (remove_chr) {
    seeksv <- seeksv %>% 
      mutate(`@left_chr` = str_replace(`@left_chr`, "chr", "")) %>% 
      mutate(right_chr = str_replace(right_chr, "chr", ""))
  }
    
  print(glue::glue("imported {nrow(seeksv)} rows, filering for integrations"))

  # get only host/virus lines
  seeksv <- seeksv %>% 
    dplyr::rename(left_chr = `@left_chr`) %>% 
    dplyr::filter(xor(stringr::str_detect(left_chr, host_chrs), 
                      stringr::str_detect(right_chr, host_chrs))) 
  
  # figure out host and virus coordinates for each event

  seeksv <- seeksv %>% 
    rowwise() %>% 
    dplyr::mutate(left_host = stringr::str_detect(left_chr, host_chrs)) %>% 
    dplyr::mutate(host_chr = ifelse(left_host, left_chr, right_chr)) %>% 
    dplyr::mutate(host_start = dplyr::case_when(
      left_host && (microhomology_length < 0) ~ left_pos + microhomology_length,
      left_host && (microhomology_length >= 0) ~ left_pos,
      !left_host && (microhomology_length < 0) ~ right_pos + microhomology_length,
      !left_host && (microhomology_length >= 0) ~ right_pos
    )) %>% 
    dplyr::mutate(host_stop= host_start + abs(microhomology_length)) %>% 
    dplyr::mutate(virus = ifelse(left_host, right_chr, left_chr)) %>% 
    dplyr::mutate(virus_start = dplyr::case_when(
      left_host && (microhomology_length < 0) ~ right_pos + microhomology_length,
      left_host && (microhomology_length >= 0) ~ right_pos,
      !left_host && (microhomology_length < 0) ~ left_pos + microhomology_length,
      !left_host && (microhomology_length >= 0) ~ left_pos
    )) %>% 
    dplyr::mutate(virus_stop = virus_start + abs(microhomology_length)) %>% 
    ungroup()
  
  print(glue::glue("imported {nrow(seeksv)} integrations"))
  
  return(seeksv)
  
}

# import seeksv results
# we need the host fai to know which chromosomes are host and which are virus
importSeeksvSRA <- function(seeksv_dir, meta, host_fai, remove_chr=FALSE) {

  seeksv <- importSeeksv(seeksv_dir, host_fai, remove_chr)
  
  return(left_join(seeksv, meta, by="Run"))
  
}

importVerse <- function(verse_dir) {
  verse_pattern <- "integration-sites.txt"
  
  coltypes <- cols(
    `Chromosome 1` = col_character(),
    `Position 1` = col_character(),
    `Strand 1` = col_character(),
    `Chromosome 2` = col_character(),
    `Position 2` = col_character(),
    `Strand 2` = col_character(),
    `#Support reads (pair+softclip)` = col_character(),
    Confidence = col_character()
  )
  
  verse <- tibble::tibble(
    f = list.files(verse_dir, verse_pattern, recursive=TRUE),
    Run = stringr::str_split(basename(dirname(f)), "\\.", simplify=TRUE)[,1],
    host = stringr::str_split(basename(dirname(f)), "\\.", simplify=TRUE)[,2],
    virus = stringr::str_split(basename(dirname(f)), "\\.", simplify=TRUE)[,3],
    data = map(f, ~read_tsv(file.path(verse_dir, .), col_types = coltypes)),
    params = "verse"
  ) %>% 
    tidyr::unnest(data)
  
  verse <- verse %>% 
    dplyr::mutate(host_chr = ifelse(`Chromosome 1` == "chrVirus", `Chromosome 2`, `Chromosome 1`)) %>% 
    dplyr::mutate(host_start = as.integer(ifelse(`Chromosome 1` == "chrVirus", `Position 2`, `Position 1`))) %>% 
    dplyr::mutate(host_stop = host_start) %>% 
    dplyr::mutate(virus = ifelse(`Chromosome 1` == "chrVirus", `Chromosome 1`, `Chromosome 2`)) %>% 
    dplyr::mutate(virus_start = ifelse(`Chromosome 1` == "chrVirus", `Position 1`, `Position 2`)) %>% 
    dplyr::mutate(virus_stop = virus_start) 
  
  return(verse)
  
}


importVerseSRA <- function(verse_dir, meta) {
  
  verse <- importVerse(verse_dir)
  
  verse <- verse %>% 
    left_join(meta, by="Run")
  
  return(verse)
  
}

# import vseq-toolkit
importVseqToolkitSRA <- function(vseq_dir, meta) {
  vseq_pattern <- "ISGenomeVector.Unclustered.csv"
  
  files <- list.files(vseq_dir, pattern = vseq_pattern, recursive=TRUE)
  
  print(glue::glue("importing {length(files)} files"))
  
  coltypes <- cols(
    Chr = col_character(),
    GenomicPosition = col_double(),
    StrandGenomic = col_character(),
    VectorName = col_character(),
    VectorPosition = col_double(),
    StrandVector = col_character(),
    ReadID = col_character(),
    Sequence = col_character(),
    SpanGenomic = col_double(),
    SpanVector = col_double(),
    OverlapFusion = col_double(),
    DistanceFusion = col_double()
  )
  
  vseq <- tibble::tibble(
    f = file.path(vseq_dir, files),
    Run = stringr::str_split(basename(dirname(f)), "\\.", simplify=TRUE)[,1],
    host = stringr::str_split(basename(dirname(f)), "\\.", simplify=TRUE)[,2],
    virus = stringr::str_split(basename(dirname(f)), "\\.", simplify=TRUE)[,3],
    data = map(f, ~read_csv(., col_types = coltypes)),
    params = "VSeq-Toolkit"
  ) %>% 
    unnest(data)
  
  return(vseq)
  
}
  

# write isling bed files for use with bedtools
# input dataframe must have 'params' and 'Run' columns
# so that we know how to separate out sites
writeBedFiles <- function(ints, dir, tool) {
  allowed_tools <- c("isling", "vifi", "polyidus", "seeksv", "verse", "vseq")
  if (!tool %in% allowed_tools) {
    stop(glue::glue("tool must be one of: {paste0(allowed_tools, collapse = ', ')}"))
  }
  
  # get unique params and runs from ints tibble
  params <- unique(ints$params)
  runs <- unique(ints$Run)
  
  for (p in params) {
    param_dir <- file.path(dir, p)
    dir.create(param_dir)
    for (r in runs) {
      
      # get data to write for this set of params and run
      if (tool == "isling") {
        to_write <- getIslingtoWrite(ints, p, r)
      } else if (tool == "polyidus") {
        to_write <- getPolyidustoWrite(ints, p, r)
      } else if (tool == "vifi") {
        to_write <- getVifitoWrite(ints, p, r)
      } else if (tool == "seeksv") {
        to_write <- getSeeksvtoWrite(ints, p, r)
      } else if (tool == "verse") {
        to_write <- getVersetoWrite(ints, p, r)
      } else if (tool == "vseq") {
        to_write <- getVseqtoWrite(ints, p, r)
      }
      
      filename <- file.path(param_dir, glue::glue("{r}.bed"))
      
      # write empty file or write data and sort
      if (nrow(to_write) == 0) {
        system(glue::glue("touch {filename}"))
      } else {
        write_tsv(to_write, filename, col_names = FALSE)
        tmp <- glue::glue("{filename}.tmp")
        system(glue::glue("sort -k1,1 -k2,2n {filename} > {tmp}"))
        system(glue::glue("mv {tmp} {filename}"))
      }
      
    }
  }
  return(params)
}

writeIslingBedFiles <- function(isling_df, dir) {

  return(writeBedFiles(isling_df, dir, "isling"))
}

getIslingtoWrite <- function(isling_df, param, run) {
  return(isling_df %>% 
           filter(params == param) %>% 
           filter(Run == run) %>% 
           select(Chr, IntStart, IntStop, ReadIDs)
  )
}

writePolyidusBedFiles <- function(polyidus, dir) {

  return(writeBedFiles(polyidus, dir, "polyidus"))
  
}

getPolyidustoWrite <- function(polyidus, param, run) {
  return(polyidus %>% 
    filter(params == param) %>% 
    filter(Run == run) %>% 
    select(chrom, chromStart, chromEnd) %>% 
    mutate(id = '.'))
}

getVifitoWrite <- function(vifi, param, run) {
  return(vifi %>% 
           filter(params == param) %>% 
           filter(Run == run) %>% 
            select(Chr, Min, Max) %>% 
            mutate(id = ".")
    )
}

writeVifiBedFiles <- function(vifi, dir) {
  
  return(writeBedFiles(vifi, dir, "vifi"))
  
}

getSeeksvtoWrite <- function(seeksv, param, run) {
  return(seeksv %>% 
           filter(params == param) %>% 
           filter(Run == run) %>% 
           select(host_chr, host_start, host_stop) %>% 
           mutate(id = ".")
         )
}

writeSeeksvBedFiles <- function(seeksv, dir) {
  return(writeBedFiles(seeksv, dir, "seeksv"))
}

getVersetoWrite <- function(verse, param, run) {
  return(
    verse %>% 
      filter(Run == run) %>% 
      filter(params == param) %>% 
      select(host_chr, host_start, host_stop) %>% 
      mutate(id = ".")
  )
}

writeVerseBedFiles <- function(verse, dir) {
  return(writeBedFiles(verse, dir, "verse"))
}

getVseqtoWrite <- function(vseq, param, run) {
  return(vseq %>% 
           filter(params == param) %>% 
           filter(Run == run) %>% 
           select(Chr, GenomicPosition) %>% 
           mutate(host_stop = GenomicPosition) %>% 
           mutate(id = '.')
         )
}
writeVseqBedFiles <- function(vseq, dir) {
  return(writeBedFiles(vseq, dir, "vseq"))
}

# sort a list of validated bed files in a folder
sortBedFiles <- function(in_dir, out_dir) {
  # list bed files in directory
  files <- list.files(in_dir, pattern = ".bed")
  print(glue::glue("found {length(files)} files for sorting"))
  
  # loop over files
  for (f in files) {
    in_file <- file.path(in_dir, f)
    file_base <- basename(tools::file_path_sans_ext(f))
    out_file <- file.path(out_dir, glue::glue("{file_base}.sorted.bed"))
    system(glue::glue("sort -k1,1 -k2,2n {in_file} > {out_file}"))
  }
}

# comapre output of various tools (in bed format - written with a function above)
# with validated or sites from authors of data
# tool output could be from runs with various parameters
# 'comparsion' means bedtools closest
compareBedFiles <- function(experimental_dir, params, validated_dir, bedtools) {
  
  validated_files <- list.files(validated_dir, pattern = ".sorted.bed")
  
  # loop over params
  for (p in params) {
    
    # get files in dir
    params_dir <- file.path(experimental_dir, p)
    if (!dir.exists(params_dir)) {
      stop(glue::glue("directory {params_dir} does not exist"))
    }
    files <- list.files(params_dir, pattern = ".bed")
    files <- files[!str_detect(files, "closest")]
    
    # loop over files in directory
    for (f in files) {
      
      # get run name
      run <- f %>% basename() %>% tools::file_path_sans_ext()
      
      # find validated file 
      a_file <- validated_files[str_detect(validated_files, run)]
      a_file <- file.path(validated_dir, a_file)
      b_file <- file.path(params_dir, f)
      
      # construct output filename
      out_file <- file.path(params_dir, glue::glue("{run}.closest.bed"))

      # do closeset
      system(glue::glue('printf "chr_ref\tstart_ref\tstop_ref\tori_ref\tint_site_ref\tChr\tIntStart\tIntStop\tid\td\n" > {out_file}'))
      cmd <- glue::glue("{bedtools} closest -d -t first -a {a_file} -b {b_file} >> {out_file}")
      system(cmd)
      
    }
  }
}

# import data from intersect
importComparison <- function(experimental_dir, params) {
  
  colnames <- cols(
    chr_ref = col_character(),
    start_ref = col_double(),
    stop_ref = col_double(),
    ori_ref = col_character(),
    int_site_ref = col_character(),
    Chr = col_character(),
    IntStart = col_double(),
    IntStop = col_double(),
    id = col_character(),
    d = col_double()
  )
  
  intersect <- tibble(
    params = params, 
    path = file.path(experimental_dir, params),
    files = map(path, ~list.files(., pattern = ".closest.bed"))
  ) %>% 
    unnest(files) %>% 
    mutate(run = str_split(basename(files), "\\.", simplify=TRUE)[,1]) %>% 
    mutate(data = map2(files, path, ~read_tsv(file.path(.y, .x), 
                                              col_types=colnames, 
                                              na = c("", "NA", '.')))) %>% 
    unnest(data)
  return(intersect)
  
}

# compare output of various tools with other tools

