#!/bin/bash

#hg19 - version recommended by Heng Li (https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use)
#wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
#gunzip human_g1k_v37.fasta.gz


#hg38 - version recommended by Heng Li (https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use)
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
#gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
#mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna hg38.fa

#hg38 refSeq
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz
#gunzip GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz
#get genes from hg38 refSeq
#awk '{ if (($0 ~ /^#/) || ($3 ~ /gene/)) { print } }' GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff > hg38_genes.gff
#get exons from hg38 refseq
#awk '{ if (($0 ~ /^#/) || ($3 ~ /exon/)) { print } }' GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff > hg38_exons.gff


#hg38 GENCODE
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gff3.gz
mv GRCh38.primary_assembly.genome.fa.gz GRCh38.fa.gz
mv gencode.v32.annotation.gff3.gz GRCh38.gencode.v32.annotation.gff3.gz


#hg19 - version used by HGT-ID
#wget ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Assembly/GRCh37-HG19_Broad_variant/Homo_sapiens_assembly19.fasta
#cat Homo_sapiens_assembly19.fasta | awk '{if($0 ~ /^>/){print ">chr"$1} else {print}}' | sed -e 's/chr>/chr/g' | sed -e 's/chrMT/chrM/g' > human.fa

#mouse GRCm38
#mkdir -p mouse
#cd mouse
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Mus_musculus/GRCm38/Primary_Assembly/assembled_chromosomes/FASTA/*.fa.gz
#cd ..
#zcat mouse/*.fa.gz | awk '{if($0 ~ /^>/){print ">chr"$5} else {print}}' | sed -e 's/chr>/chr/g' | sed -e 's/,//g' > mouse.fa
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz
#gunzip GCF_000001635.26_GRCm38.p6_genomic.fna.gz
#mv GCF_000001635.26_GRCm38.p6_genomic.fna mouse.fa
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.gff.gz
#gunzip *gz
#mv GCF_000001635.26_GRCm38.p6_genomic.fna mouse.fa

#mouse genome from GENCODE
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gff3.gz
mv GRCm38.primary_assembly.genome.fa.gz GRCm38.fa.gz
mv gencode.vM23.annotation.gff3.gz GRCm38.gencode.vM23.annotation.gff3.gz
gunzip *.gz


#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz 
#gunzip GCF_000001635.26_GRCm38.p6_genomic.fna.gz
#cat GCF_000001635.26_GRCm38.p6_genomic.fna | awk '{if($0 ~ /^>/){print ">chr"$1} else {print}}'  | sed -e 's/chr>/chr/g' > mouse.fa

## macaque macaca fasicularis
mkdir -p macaque
cd macaque
wget ftp://ftp.ensembl.org/pub/release-97/fasta/macaca_fascicularis/dna/*dna.chromosome*.gz

#zcat *.fa.gz | perl -pe 's/^>(\w+).+$/>chr$1/' | perl -pe 's/^>chrMT/>chrM/' > ../macaque.fa
zcat *.fa.gz > ../macaque.fa

cd ..

wget ftp://ftp.ensembl.org/pub/release-97/gff3/macaca_fascicularis/Macaca_fascicularis.Macaca_fascicularis_5.0.97.gff3.gz
gunzip Macaca_fascicularis.Macaca_fascicularis_5.0.97.gff3.gz
#sed s/^/chr/g Macaca_fascicularis.Macaca_fascicularis_5.0.97.gff3 | sed 's/^chr##/##/g' > Macaca_fascicularis_5.0.97.gff
mv Macaca_fascicularis.Macaca_fascicularis_5.0.97.gff3 Macaca_fascicularis_5.0.97.gff

## to compare against sunando's results, use MacFas5 from genbank
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Macaca_fascicularis/latest_assembly_versions/GCA_000364345.1_Macaca_fascicularis_5.0/GCA_000364345.1_Macaca_fascicularis_5.0_genomic.fna.gz
gunzip GCA_000364345.1_Macaca_fascicularis_5.0_genomic.fna.gz
mv GCA_000364345.1_Macaca_fascicularis_5.0_genomic.fna macFas5_genbank.fa
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Macaca_fascicularis/latest_assembly_versions/GCA_000364345.1_Macaca_fascicularis_5.0/GCA_000364345.1_Macaca_fascicularis_5.0_genomic.gff.gz
gunzip GCA_000364345.1_Macaca_fascicularis_5.0_genomic.gff.gz


#### to generate a bed file of gaps (streches of N's) in a fasta file:
# https://www.biostars.org/p/133742/#377084

perl -ne 'chomp;if( />(.*)/){$head = $1; $i=0; next};@a=split("",$_); foreach(@a){$i++; if($_ eq "N" && $s ==0 ){$z=$i-1; print "$head\t$z"; $s =1}elsif($s==1 && $_ ne "N"){$j=$i-1;print "\t$j\n";$s=0}}' ref.fa > ref_gaps.bed

#depending on the format of the fasta file, might need to remove some columns
awk '{print $1"\t"$5"\t"$6}' macFas5_gaps.bed > macFas5_gaps_.bed



