#!/bin/bash

#hg19 - version recommended by Heng Li (https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use)
#wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
#gunzip human_g1k_v37.fasta.gz


#hg38 - version recommended by Heng Li (https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna hg38.fa

#hg38 refSeq
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz


#hg19 - version used by HGT-ID
wget ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Assembly/GRCh37-HG19_Broad_variant/Homo_sapiens_assembly19.fasta
cat Homo_sapiens_assembly19.fasta | awk '{if($0 ~ /^>/){print ">chr"$1} else {print}}' | sed -e 's/chr>/chr/g' | sed -e 's/chrMT/chrM/g' > human.fa

#mouse GRCm38
#mkdir -p mouse
#cd mouse
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Mus_musculus/GRCm38/Primary_Assembly/assembled_chromosomes/FASTA/*.fa.gz
#cd ..
#zcat mouse/*.fa.gz | awk '{if($0 ~ /^>/){print ">chr"$5} else {print}}' | sed -e 's/chr>/chr/g' | sed -e 's/,//g' > mouse.fa
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.gff.gz
gunzip *gz
mv GCF_000001635.26_GRCm38.p6_genomic.fna mouse.fa


#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz 
#gunzip GCF_000001635.26_GRCm38.p6_genomic.fna.gz
#cat GCF_000001635.26_GRCm38.p6_genomic.fna | awk '{if($0 ~ /^>/){print ">chr"$1} else {print}}'  | sed -e 's/chr>/chr/g' > mouse.fa



