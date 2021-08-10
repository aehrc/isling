#!/bin/bash
set -euo pipefail

cd data/references/

#hg37
if [ ! -e human_g1k_v37.fasta ] ; then
	echo "downloading human reference (GRCh37/hg19)"
	wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/
	gunzip human_g1k_v37.fasta.gz
fi

if [ ! -e human_g1k_v37.fasta.fai ] ; then
	samtools faidx human_g1k_v37.fasta
fi

#mm10
if [ ! -e mm10_no_alt_analysis_set_ENCODE.fasta ] ; then
	echo "downloading mm10 reference"
	wget https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz
	gunzip mm10_no_alt_analysis_set_ENCODE.fasta.gz
fi

if [ ! -e mm10_no_alt_analysis_set_ENCODE.fasta.fai ] ; then
	samtools faidx mm10_no_alt_analysis_set_ENCODE.fasta
fi
