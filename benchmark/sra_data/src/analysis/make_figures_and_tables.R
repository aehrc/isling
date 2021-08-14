
#### imports ####


library(tidyverse)
library(scales)
source("support_funcs.R")

#### paths and setup ####

bedtools <- "/apps/bedtools/2.29.2/bin/bedtools"   

nelson_metadata_dir <- "../../data/metadata/PRJNA485509/"
nelson_sra_run_table <- file.path(nelson_metadata_dir, "SraRunTable.txt")

nelson_isling_dir <- "../../out/Nelson_2019_PRJNA485509/isling/"

nelson_polyidus_dir <- "../../out/Nelson_2019_PRJNA485509/polyidus/"

nelson_seeksv_dir <- "../../out/Nelson_2019_PRJNA485509/seeksv/ints"
nelson_host_fai <- "../../data/references/mm10_no_alt_analysis_set_ENCODE.fasta.fai"

nelson_vseq_dirs <- c("../../out/Nelson_2019_PRJNA485509/vseq_toolkit/",
                      "../../out/Nelson_2019_PRJNA485509_vseq_novecvecfusion/vseq_toolkit/")

nelson_vseq_suffix <- "ISGenomeVector.UniqueGenome.csv"

# for plotting
colors <- list(
  "isling" = "#F8766D",
  "Polyidus" = "#A3A500",
  "seeksv" = "#00BF7D",
  "ViFi" = "#00B0F6",
  "VSeq-Toolkit" = "#E76BF3"
)

dir.create("../../tables")
dir.create("../../figures")



#### nelson data - metadata ####

print("starting analysis of Nelson data")
nelson_meta <- read_tsv(nelson_sra_run_table)

#### nelson data - isling ####

print("starting analysis of Nelson results for isling")

# import data
nelson_isling <- importIslingHostSRA(nelson_isling_dir, nelson_meta)

nelson_isling

# compare results to data from original paper

tmpdir <- "tmp"
system(glue::glue("rm -rf {tmpdir}"))
dir.create(tmpdir)

nelson_isling_write <- nelson_isling %>% 
  mutate(Run = "AAV_integrations")

params <- writeIslingBedFiles(nelson_isling_write, tmpdir)

validated <- "tmp/validated"
dir.create(validated)
sortBedFiles(nelson_metadata_dir , validated)

compareBedFiles(tmpdir, params, validated, bedtools)

#params <- c("all", "filter1", "filter2", "filter3")
nelson_isling_intersects <- importComparison(tmpdir, params)

#### nelson data- polyidus ####

print("starting analysis of Nelson results for polyidus")

# import
nelson_polyidus <- importPolyidusSRA(nelson_polyidus_dir, nelson_meta)

# write to file for bedtools
nelson_poly_write <- nelson_polyidus %>% 
  mutate(Run = "AAV_integrations")
params <- writePolyidusBedFiles(nelson_poly_write, tmpdir)

# compare polyidus sites against published sites
compareBedFiles(tmpdir, params, validated, bedtools)

# import results
nelson_poly_intersects <- importComparison(tmpdir, params)

#### nelson data - seeksv ####

print("starting analysis of Nelson results for seeksv")

# import seeksv
nelson_seeksv <- importSeeksvSRA(nelson_seeksv_dir, nelson_meta, nelson_host_fai) %>% 
  filter(str_detect(left_chr, "paragRNA1") | str_detect(right_chr, "paragRNA1"))

# compare against published sites
nelson_seeksv_write <- nelson_seeksv %>% 
  mutate(Run = "AAV_integrations")

params <- writeSeeksvBedFiles(nelson_seeksv_write, tmpdir)

compareBedFiles(tmpdir, params, validated, bedtools)

nelson_seeksv_intersects <- importComparison(tmpdir, params)

#### nelson data - vseq-toolkit ####

print("starting analysis of Nelson results for vseq-toolkit")

# import two different analysis parameters
nelson_vseq_dir <- nelson_vseq_dirs[1]
nelson_vseq <- importVseqToolkitSRA(nelson_vseq_dir, nelson_meta) %>% 
  mutate(params = "VSeq_vec_vec_fusion")

print(nelson_vseq_dirs[2])
nelson_vseq_dir <- nelson_vseq_dirs[2]
nelson_vseq2 <- importVseqToolkitSRA(nelson_vseq_dir, nelson_meta) %>% 
  mutate(params = "VSeq_no_vec_vec_fusion")

nelson_vseq <- bind_rows(nelson_vseq, nelson_vseq2) %>% 
  filter(str_detect(VectorName, "paragRNA1"))

# compoare with published
nelson_vseq_write <- nelson_vseq %>% 
  mutate(Run = "AAV_integrations")

params <- writeVseqBedFiles(nelson_vseq_write, tmpdir)

compareBedFiles(tmpdir, params, validated, bedtools)

nelson_vseq_intersects <- importComparison(tmpdir, params)

#### nelson - all tools combined ####

print("combining all results for Nelson data")

nelson_to_join <- list(nelson_isling_intersects, nelson_poly_intersects, nelson_vseq_intersects, nelson_seeksv_intersects)
nelson_combined <- tibble::tibble()

for (result in nelson_to_join) {
  nelson_combined <- bind_rows(nelson_combined, result)
}

nelson_results <- nelson_combined %>% 
  distinct(params, run, chr_ref, start_ref, stop_ref, int_site_ref, d) %>% 
  pivot_wider(names_from = "params", values_from = "d") %>% 
  mutate(int_site_ref = as.numeric(int_site_ref)) %>% 
  arrange(int_site_ref)

nelson_results <- nelson_results %>% 
  mutate(dataset = "Nelson")

dir.create("../../tables")
dir.create("../../figures")

nelson_results %>% 
  select(-run) %>% 
  rename(Chr = chr_ref) %>% 
  rename(Pos = start_ref) %>% 
  rename(isling = filter3) %>% 
  rename(Polyidus=polyidus) %>% 
  rename(VSeq_Toolkit = VSeq_no_vec_vec_fusion) %>% 
  rename(seeksv = seeksv) %>% 
  select(Chr, Pos, isling, Polyidus, VSeq_Toolkit, seeksv) %>% 
  write_tsv("../../tables/supp-table-7.tsv")


# make plot
nelson_results_filt <- nelson_results %>% 
  pivot_longer(all:seeksv, names_to="tool", values_to="distance") %>% 
  mutate( tool = case_when(
    tool == "filter3" ~ "isling",
    tool == "polyidus" ~ "Polyidus",
    #    tool == "VSeq_vec_vec_fusion" ~ "VSeq-Toolkit (removing Vec-vec fusion)",
    tool == "VSeq_no_vec_vec_fusion" ~ "VSeq-Toolkit",
    tool == "seeksv" ~ "seeksv"
  )) %>% 
  filter(!is.na(tool))

nelson_freqpoly <- nelson_results_filt %>% 
  mutate(distance = distance + 0.1) %>% 
  ggplot(aes(x = distance, color=tool)) +
  geom_freqpoly(bins=40) +
  scale_x_log10() +
  scale_color_manual(values=c("#A3A500", "#E76BF3", "#F8766D", "#00BF7D"))

nelson_freqpoly

# threshold to calculate found/missed integrations
nelson_n_plot <- nelson_results_filt %>% 
  mutate(found_5 = case_when(
    distance < 0 ~ FALSE,
    distance <= 5 ~ TRUE,
    distance > 5 ~ FALSE
  )) %>% 
  mutate(found_100 = case_when(
    distance < 0 ~ FALSE,
    distance <= 100 ~ TRUE,
    distance > 100 ~ FALSE
  )) %>%   
  pivot_longer(found_5:found_100,
               names_to="threshold", names_prefix="found_", values_to="found") %>% 
  group_by(tool, run, threshold) %>% 
  summarise(not_found = sum(!found),
            found = sum(found)) %>% 
  ungroup() %>% 
  pivot_longer(not_found:found, names_to="type", values_to="number") %>%  
  mutate(tool = forcats::fct_relevel(tool, "isling", "Polyidus", "seeksv", "VSeq-Toolkit")) %>% 
  mutate(type = forcats::fct_relevel(type, "not_found", "found")) %>% 
  mutate(threshold = forcats::fct_relevel(as.character(threshold),"5", "100")) %>% 
  ggplot(aes(x = tool, y = number, fill=type)) +
  geom_col() +
  scale_fill_manual(values = c("grey", "black")) +
  facet_wrap(~threshold)

nelson_n_plot


#### Sung data - paths ####

print("starting analysis of Sung results")

sung_metadata_dir <- "../../data/metadata/PRJEB2869"
sung_sra_run_table <- file.path(sung_metadata_dir, "SraRunTable.txt")

wang_accs <- file.path(sung_metadata_dir, "SRR_Acc_List.txt")
sung_bed_dir <- file.path(sung_metadata_dir, "validated")

sung_isling_dir <- "../../out/Sung_2012_PRJEB2869_hg19/isling/"

sung_polyidus_dir <- "../../out/Sung_2012_PRJEB2869_hg19/polyidus"

sung_vifi_dir <- "../../out/Sung_2012_PRJEB2869_hg19/vifi"

sung_seeksv_dir <- "../../out/Sung_2012_PRJEB2869_hg19/seeksv/ints"
sung_host_fai <- "../../data/references/human_g1k_v37.fasta.fai"

#### Sung - metadata ####

sung_meta <- read_csv(sung_sra_run_table)

# only have data for these accessions - only the ones with integrations
accs <- read_tsv(wang_accs, col_names=FALSE)

sung_meta <- sung_meta %>% 
  filter(Run %in% accs$X1)

#### Sung - isling ####

print("starting analysis of Sung results for isling")

# import data
sung_isling <- importIslingHostSRA(sung_isling_dir, sung_meta)

# compoare with published
tmpdir <- "tmp"
system(glue::glue("rm -rf {tmpdir}"))
dir.create(tmpdir)

sung_isling_write <- sung_isling %>% 
  mutate(Run = Alias)

params <- writeIslingBedFiles(sung_isling_write, tmpdir)

validated <- "tmp/validated"
dir.create(validated)
sortBedFiles(sung_bed_dir , validated)

compareBedFiles(tmpdir, params, validated, bedtools)

sung_isling_intersects <- importComparison(tmpdir, params)


# include all integrations (not just unambiguous ones) and compare
sung_isling_all <- importIslingAllSRA(sung_isling_dir, sung_meta)

sung_isling_all_write <- sung_isling_all %>% 
  mutate(Run = Alias)

params <- writeIslingBedFiles(sung_isling_all_write, tmpdir)

compareBedFiles(tmpdir, params, validated, bedtools)

sung_isling_intersects_all <- importComparison(tmpdir, params)

#### Sung - polyidus ####

print("starting analysis of Sung results for polyidus")

sung_polyidus <- importPolyidus(sung_polyidus_dir)

# compare against published sites
params <- writePolyidusBedFiles(sung_polyidus, tmpdir)

compareBedFiles(tmpdir, params, validated, bedtools)

sung_poly_intersects <- importComparison(tmpdir, params)

#### sung - ViFi ####

print("starting analysis of Sung results for vifi")

sung_vifi <- importVifi(sung_vifi_dir)

# compare against published sites
params <- writeVifiBedFiles(sung_vifi, tmpdir)

compareBedFiles(tmpdir, params, validated, bedtools)

sung_vifi_intersects <- importComparison(tmpdir, params)

#### sung - seeksv ####

print("starting analysis of Sung results for seeksv")

sung_seeksv <- importSeeksv(sung_seeksv_dir,sung_host_fai)

params <- writeSeeksvBedFiles(sung_seeksv, tmpdir)

compareBedFiles(tmpdir, params, validated, bedtools)

sung_seeksv_intersects <- importComparison(tmpdir, params)

#### Sung - combine all tools ####

print("combining all results for Sung data")

sung_to_join <- list(sung_isling_intersects, sung_isling_intersects_all, sung_poly_intersects, sung_vifi_intersects, sung_seeksv_intersects)
sung_combined <- tibble::tibble()

for (result in sung_to_join) {
  sung_combined <- bind_rows(sung_combined, result)
}

sung_results <- sung_combined %>% 
  distinct(params, run, chr_ref, start_ref, stop_ref, int_site_ref, d) %>% 
  pivot_wider(names_from = "params", values_from = "d") %>% 
  mutate(gene = str_split(int_site_ref, "_", simplify=TRUE)[,3])

sung_results %>% 
  rename(Sample = run) %>% 
  rename(Chr = chr_ref) %>% 
  rename(Pos = start_ref) %>% 
  rename(Gene = gene) %>% 
  rename(isling = all) %>% 
  rename(`isling (including ambiguous location)` = `all_unique-plus-ambig`) %>% 
  rename(Polyidus = polyidus) %>% 
  rename(ViFi = vifi) %>% 
  select(Sample, Chr, Pos, Gene, isling, `isling (including ambiguous location)`, Polyidus, ViFi) %>% 
  write_tsv("../../tables/supp-table-5.tsv")

sung_results_filt <- sung_results %>% 
  pivot_longer(all:seeksv, names_to="tool", values_to="distance") %>% 
  mutate( tool = case_when(
    tool == "all" ~ "isling",
    tool == "polyidus" ~ "Polyidus",
    tool == "seeksv" ~ "seeksv",
    tool == "vifi" ~ "ViFi"
  )) %>% 
  filter(!is.na(tool))

# count integrations within thrshold distance of published site
sung_results_filt <- sung_results %>% 
  pivot_longer(all:seeksv, names_to="tool", values_to="distance") %>% 
  mutate( tool = case_when(
    tool == "all" ~ "isling",
    tool == "polyidus" ~ "Polyidus",
    tool == "seeksv" ~ "seeksv",
    tool == "vifi" ~ "ViFi"
  )) %>% 
  filter(!is.na(tool))
dist_threshold <- 5

sung_n_plot <- sung_results_filt %>% 
  mutate(found = case_when(
    distance < 0 ~ FALSE,
    distance <= dist_threshold ~ TRUE,
    distance > dist_threshold ~ FALSE,
    is.na(distance) ~ FALSE,
  )) %>% 
  group_by(tool) %>% 
  summarise(not_found = sum(!found),
            found = sum(found)) %>% 
  ungroup() %>% 
  pivot_longer(not_found:found, names_to="type", values_to="number") %>%  
  mutate(tool = forcats::fct_relevel(tool, "isling", "Polyidus", "ViFi", "seeksv" )) %>% 
  mutate(type = forcats::fct_relevel(type, "not_found", "found")) %>% 
  ggplot(aes(x = tool, y = number, fill=type)) +
  geom_col() +
  scale_fill_manual(values = c("grey", "black"))

sung_n_plot 

#### Lau - paths ####

print("starting analysis of Lau data")

lau_metadata_dir <- "../../data/metadata/SRP023539"
lau_sra_run_table <- file.path(lau_metadata_dir, "SraRunTable.txt")

lau_isling_dir <- "../../out/Lau_2014_SRP023539_hg19/isling/"

lau_polyidus_dir <- "../../out/Lau_2014_SRP023539_hg19/polyidus/"

lau_vifi_dir <- "../../out/Lau_2014_SRP023539_hg19/vifi/"

lau_host_fai <- "../../data/references/human_g1k_v37.fasta.fai"
lau_virus_fai <- "../../data/references/SRP023539_HBV.fa.fai"
lau_seeksv_dir <- "../../out/Lau_2014_SRP023539_hg19/seeksv/ints"

lau_vseq_dirs <- c("../../out/Lau_2014_SRP023539_hg19/vseq_toolkit/")
lau_vseq_suffix <- "ISGenomeVector.UniqueGenome.csv"

#### lau - metadata ####

lau_meta <- read_csv(lau_sra_run_table)

#### lau - isling ####

print("starting analysis of Lau results for isling")

lau_isling <- importIslingHostSRA(lau_isling_dir, lau_meta)

# write bed files from each set of parameters and each run
tmpdir <- "tmp"
system(glue::glue("rm -r {tmpdir}"))
dir.create(tmpdir)
params <- writeIslingBedFiles(lau_isling, tmpdir)

validated <- "tmp/validated"
dir.create(validated)
sortBedFiles(lau_metadata_dir, validated)

compareBedFiles(tmpdir, params, validated, bedtools)

lau_isling_intersects <- importComparison(tmpdir, params)

# include all integrations (not just unambiguous ones) and compare
lau_isling_all <- importIslingAllSRA(lau_isling_dir, lau_meta)

lau_isling_all_write <- lau_isling_all

params <- writeIslingBedFiles(lau_isling_all_write, tmpdir)

compareBedFiles(tmpdir, params, validated, bedtools)

lau_isling_intersects_all <- importComparison(tmpdir, params)

#### lau - polyidus ####

print("starting analysis of Lau results for Polyidus")

lau_polyidus <- importPolyidusSRA(lau_polyidus_dir, lau_meta)

params <- writePolyidusBedFiles(lau_polyidus, tmpdir)

compareBedFiles(tmpdir, params, validated, bedtools)

lau_poly_intersects <- importComparison(tmpdir, params)

#### lau - vifi ####

print("starting analysis of Lau results for ViFi")

lau_vifi <- importVifiSRA(lau_vifi_dir, lau_meta)

params <- writeVifiBedFiles(lau_vifi, tmpdir)

compareBedFiles(tmpdir, params, validated, bedtools)

lau_vifi_intersects <- importComparison(tmpdir, params)

#### lau - seeksv ####

print("starting analysis of Lau results for seeksv")

lau_seeksv <- importSeeksvSRA(lau_seeksv_dir, lau_meta, lau_host_fai)

params <- writeSeeksvBedFiles(lau_seeksv, tmpdir)

compareBedFiles(tmpdir, params, validated, bedtools)

lau_seeksv_intersects <- importComparison(tmpdir, params)


#### lau - combined all tools ####

print("combining all results for Lau data")

lau_to_join <- list(lau_isling_intersects, lau_isling_intersects_all, lau_poly_intersects, lau_vifi_intersects,  lau_seeksv_intersects)
lau_combined <- tibble::tibble()

for (result in lau_to_join) {
  lau_combined <- bind_rows(lau_combined, result)
}

lau_results <- lau_combined %>% 
  distinct(params, run, chr_ref, start_ref, stop_ref, int_site_ref, d) %>% 
  pivot_wider(names_from = "params", values_from = "d") %>% 
  mutate(cell_line = str_split(int_site_ref, "_", simplify=TRUE)[,1])

lau_results %>% 
  rename(`Cell Line` = cell_line) %>% 
  rename(`SRA Run` = run) %>% 
  rename(Chr = chr_ref) %>% 
  rename(Start = start_ref) %>% 
  rename(Stop = stop_ref) %>% 
  rename(isling = filter1) %>% 
  rename(`isling (including ambiguous)` = `filter1_unique-plus-ambig`) %>% 
  rename(Polyidus = polyidus) %>% 
  rename(ViFi = vifi) %>% 
  select(`Cell Line`, `SRA Run`, `Chr`, `Start`, `Stop`,  isling, `isling (including ambiguous)`, Polyidus, ViFi, seeksv) %>% 
  write_tsv("../../tables/supp-table-6.tsv")

lau_results_filt <- lau_results %>% 
  pivot_longer(filter1:seeksv, names_to="tool", values_to="distance") %>% 
  mutate( tool = case_when(
    tool == "filter1" ~ "isling",
    tool == "polyidus" ~ "Polyidus",
    tool == "seeksv" ~ "seeksv",
    tool == "vifi" ~ "ViFi"
  )) %>% 
  filter(!is.na(tool))

dist_threshold <- 5

lau_n_plot <- lau_results_filt %>% 
  mutate(found = case_when(
    distance < 0 ~ FALSE,
    distance <= dist_threshold ~ TRUE,
    distance > dist_threshold ~ FALSE,
    is.na(distance) ~ FALSE,
  )) %>% 
  group_by(tool) %>% 
  summarise(not_found = sum(!found),
            found = sum(found)) %>% 
  ungroup() %>% 
  pivot_longer(not_found:found, names_to="type", values_to="number") %>%  
  mutate(tool = forcats::fct_relevel(tool, "isling", "Polyidus", "ViFi", "seeksv" )) %>% 
  mutate(type = forcats::fct_relevel(type, "not_found", "found")) %>% 
  ggplot(aes(x = tool, y = number, fill=type)) +
  geom_col() +
  scale_fill_manual(values = c("grey", "black"))

lau_n_plot 

#### combined figure ####

print("generating combined figure")

type_legend <- cowplot::get_legend(lau_n_plot + theme(legend.position = "bottom"))
color_legend <- cowplot::get_legend(nelson_freqpoly + theme(legend.position = "bottom"))

sung_n_plot <- sung_n_plot + 
  theme_classic() +
  ylab("count") +
  theme(axis.text.x = element_text(angle = 45,  hjust=1),
        legend.position = "none",
        axis.title.x = element_blank()) +
  ggtitle("HBV, WGS")
nelson_n_plot <- nelson_n_plot  + 
  theme_classic() +
  theme(legend.position = "none") +
  ylab("count") +
  theme(axis.text.x = element_text(angle = 45,  hjust=1),
        axis.title.x = element_blank(),
        strip.background = element_blank()) +
  ggtitle("rAAV vector")
lau_n_plot <- lau_n_plot + 
  theme_classic() +
  ylab("count") +
  theme(axis.text.x = element_text(angle = 45,  hjust=1),
        legend.position = "none",
        axis.title.x = element_blank()) +
  ggtitle("HBV, RNA-seq")
nelson_freqpoly <- nelson_freqpoly +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,  hjust=1),
        legend.position = "none")  +
  ggtitle("rAAV vector") +
  geom_vline(xintercept=5) +
  geom_vline(xintercept=100, linetype="dotted")



plots <- cowplot::plot_grid(
  sung_n_plot,
  lau_n_plot,
  nelson_n_plot,
  nelson_freqpoly,
  labels="AUTO"
)
legends <- cowplot::plot_grid(type_legend,color_legend)

real_data <- cowplot::plot_grid(
  plots, legends,
  nrow = 2,
  rel_heights = c(1, 0.05)
)

print(real_data)

cowplot::save_plot("../../figures/figure4.pdf", real_data)


#### combined table of number of integrations found ####


thresh <- 5
lau_found_5 <- lau_results %>% 
  pivot_longer(all:seeksv, names_to = "condition", values_to = "distance") %>% 
  mutate(found = case_when(
    distance > thresh ~ FALSE,
    distance < 0 ~ FALSE,
    TRUE ~ TRUE
  )) %>% 
  group_by(condition) %>% 
  summarise(Accession = "PRJNA298941",
            Host = "Human (hg19)",
            Virus = "Hepatitis B",
            `Library type` = "RNA-seq",
            `Distance threshold` = thresh, 
            `Total sites` = n(),
            n_found = sum(found)
            ) %>% 
  pivot_wider(names_from = condition, values_from = n_found) %>% 
  rename(isling = filter1) %>% 
  rename(`isling (including ambiguous locations)` = `filter1_unique-plus-ambig`) %>% 
  rename(Polyidus = polyidus) %>% 
  rename(ViFi = vifi) %>% 
  select(Accession, Host, Virus, `Library type`, `Distance threshold`, 
         `Total sites`, isling, `isling (including ambiguous locations)`,
         Polyidus, seeksv, ViFi)
  
lau_found_5

thresh <- 8
lau_found_8 <- lau_results %>% 
  pivot_longer(all:seeksv, names_to = "condition", values_to = "distance") %>% 
  mutate(found = case_when(
    distance > thresh ~ FALSE,
    distance < 0 ~ FALSE,
    TRUE ~ TRUE
  )) %>% 
  group_by(condition) %>% 
  summarise(Accession = "PRJNA298941",
            Host = "Human (hg19)",
            Virus = "Hepatitis B",
            `Library type` = "RNA-seq",
            `Distance threshold` = thresh, 
            `Total sites` = n(),
            n_found = sum(found)
  ) %>% 
  pivot_wider(names_from = condition, values_from = n_found) %>% 
  rename(isling = filter1) %>% 
  rename(`isling (including ambiguous locations)` = `filter1_unique-plus-ambig`) %>% 
  rename(Polyidus = polyidus) %>% 
  rename(ViFi = vifi) %>% 
  select(Accession, Host, Virus, `Library type`, `Distance threshold`, 
         `Total sites`, isling, `isling (including ambiguous locations)`,
         Polyidus, seeksv, ViFi)

lau_found_8

thresh <- 5
sung_found_5 <- sung_results %>% 
  pivot_longer(all:seeksv, names_to = "condition", values_to = "distance") %>% 
  mutate(found = case_when(
    distance > thresh ~ FALSE,
    distance < 0 ~ FALSE,
    TRUE ~ TRUE
  )) %>% 
  group_by(condition) %>% 
  summarise(Accession = "PRJEB2869",
            Host = "Human (hg19)",
            Virus = "Hepatitis B",
            `Library type` = "WGS",
            `Distance threshold` = thresh, 
            `Total sites` = n(),
            n_found = sum(found)
  ) %>% 
  pivot_wider(names_from = condition, values_from = n_found)  %>% 
  rename(isling = all) %>% 
  rename(`isling (including ambiguous locations)` = `all_unique-plus-ambig`) %>% 
  rename(Polyidus = polyidus) %>% 
  rename(ViFi = vifi) %>% 
  select(Accession, Host, Virus, `Library type`, `Distance threshold`, 
         `Total sites`, isling, `isling (including ambiguous locations)`,
         Polyidus, seeksv, ViFi)
  
sung_found_5

thresh <- 5
nelson_found_5 <- nelson_results  %>% 
  pivot_longer(all:seeksv, names_to = "condition", values_to = "distance") %>% 
  mutate(found = case_when(
    distance > thresh ~ FALSE,
    distance < 0 ~ FALSE,
    TRUE ~ TRUE
  )) %>% 
  group_by(condition) %>% 
  summarise(Accession = "PRJNA485509",
            Host = "Mouse (mm10)",
            Virus = "Dual AAV/CRISPR therapy",
            `Library type` = "Amplicon / Nextera",
            `Distance threshold` = thresh, 
            `Total sites` = n(),
            n_found = sum(found)
  ) %>% 
  pivot_wider(names_from = condition, values_from = n_found) %>% 
  rename(isling = filter1) %>% 
  rename(Polyidus = polyidus) %>% 
  rename(`VSeq-Toolkit`=VSeq_no_vec_vec_fusion) %>% 
  select(Accession, Host, Virus, `Library type`, `Distance threshold`, 
         `Total sites`, isling,
         Polyidus, seeksv)


thresh <- 100
nelson_found_100 <- nelson_results  %>% 
  pivot_longer(all:seeksv, names_to = "condition", values_to = "distance") %>% 
  mutate(found = case_when(
    distance > thresh ~ FALSE,
    distance < 0 ~ FALSE,
    TRUE ~ TRUE
  )) %>% 
  group_by(condition) %>% 
  summarise(Accession = "PRJNA485509",
            Host = "Mouse (mm10)",
            Virus = "Dual AAV/CRISPR therapy",
            `Library type` = "Amplicon / Nextera",
            `Distance threshold` = thresh, 
            `Total sites` = n(),
            n_found = sum(found)
  ) %>% 
  pivot_wider(names_from = condition, values_from = n_found) %>% 
  rename(isling = filter1) %>% 
  rename(Polyidus = polyidus) %>% 
  rename(`VSeq-Toolkit`=VSeq_no_vec_vec_fusion) %>% 
  select(Accession, Host, Virus, `Library type`, `Distance threshold`, 
         `Total sites`, isling,
         Polyidus, seeksv)

summary <- bind_rows(
  sung_found_5,
  lau_found_5,
  lau_found_8,
  nelson_found_5,
  nelson_found_100
)

summary %>% 
  write_tsv("../../tables/supp-table-4.tsv")
