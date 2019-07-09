#!/bin/bash

pwd; hostname; date

#need seqprep via conda
source ~/.bashrc
conda activate integration

#get info from inputs
DATADIR=$(echo $1 | cut -f1 -d,)
DATA="${DATADIR}/reads"
SAMPLE=$(echo $1 | cut -f2 -d,)
HOST=$(echo $1 | cut -f3 -d,)
VIRUS=$(echo $1 | cut -f4 -d,)
PROJ=$2
OUTPATH="${DATADIR}/merged_reads"

#get paths to reads
READ1=${DATA}/${SAMPLE}_R1.fastq.gz
READ2=${DATA}/${SAMPLE}_R2.fastq.gz

#generate paths for output
OUT1=${OUTPATH}/${SAMPLE}.seqPrep_processed.R1.fastq.gz
OUT2=${OUTPATH}/${SAMPLE}.seqPrep_processed.R2.fastq.gz
OUT3=${OUTPATH}/${SAMPLE}.seqPrep_processed.merged.fastq.gz
OUT=${OUTPATH}/${SAMPLE}.seqPrep_processed.fastq

SeqPrep -f ${READ1} -r ${READ2} -1 ${OUT1} -2 ${OUT2} -s ${OUT3}

gzip -cdk ${OUTPATH}/${SAMPLE}* > ${OUT}


date

