
#run from directory from below 'post'
# ie Rscript post/test_postprocess.R

call <- "Rscript post/postprocess.R ../out/FRG_OTC/ints/FRG203_C73YB_CGAACTTA.GRCh38.pAAV2-OTC.integrations.txt"
dedup <- "dedup"
filter <- "filter"
mask_exclude <- c("mask-exclude", "../data/references/GRCh38/hg38_homologoustoOTCvec.bed")
nearest_gtf <- c("nearest-gtf", "../data/references/GRCh38/GENCODE/GRChg38.gencode.v33.annotation.genes.sorted.gtf")
nearest_bed <- c("nearest-bed", "../data/references/GRCh38/ENCODE/histone_ChIP-seq/ENCFF033CCP.sorted.bed")
RNA_seq <- c("RNA-seq-gtf", "../data/references/GRCh38/GENCODE/GRChg38.gencode.v33.annotation.genes.sorted.gtf", "TPM", "../data/references/GRCh38/ENCODE/RNA-seq/ENCFF143IYG.tsv")
RNA_seq2 <- c("RNA-seq-gtf", "../data/references/GRCh38/GENCODE/GRChg38.gencode.v33.annotation.genes.sorted.gtf", "TPM", "../data/references/GRCh38/ENCODE/RNA-seq/ENCFF973HRL.tsv")

system_call <- paste(call, 
                      dedup,
                      filter, 
                      paste0(mask_exclude, collapse = " "),
                      paste0(nearest_gtf, collapse = " "),
                      paste0(nearest_bed, collapse = " "),
                      paste0(RNA_seq, collapse = " "),
                      paste0(RNA_seq2, collapse = " "),
                      collapse = " ")

print(system_call)

system(system_call)

