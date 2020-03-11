
#setwd("post/")

call <- "Rscript postprocess.R ../../out/FRG_OTC/ints/FRG203_C73YB_CGAACTTA_L001.GRCh38.pAAV2-OTC.integrations.txt"
dedup <- "dedup"
filter <- "filter"
mask_exclude <- c("mask-exclude", "../../data/references/GRCh38/hg38_homologoustoOTCvec.bed ")

system_call <- paste(call, 
                      dedup,
                      filter, 
                      paste0(mask_exclude, collapse = " "),
                      collapse = " ")

system(system_call)

