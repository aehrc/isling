#!/bin/bash


# This script aligns merged reads to human and viral genomes
# Pass in (path to data, sample name, host reference, virus reference) as first argument, project directory as second argument

pwd ; hostname ; date

#get info from inputs
DATADIR=$(echo $1 | cut -f1 -d,)
DATA="${DATADIR}/merged_reads"
SAMPLE=$(echo $1 | cut -f2 -d,)
HOSTORG=$(echo $1 | cut -f3 -d,)
VIRUS=$(echo $1 | cut -f4 -d,)
PROJ=$2
BWAPATH=$3

source ~/.bashrc
conda activate integration

#input paths

#location of reads 
READS=${DATA}/${SAMPLE}.seqPrep_processed.fastq
READ1=${DATA}/${SAMPLE}.seqPrep_processed.R1.fastq.gz
READ2=${DATA}/${SAMPLE}.seqPrep_processed.R2.fastq.gz

#output paths
HOSTPATH=${DATADIR}/aligned
VIRALPATH=${DATADIR}/aligned


#index locations
HOSTINDEX=${PROJ}/references/${HOSTORG}
VIRALINDEX=${PROJ}/references/${VIRUS}


#output files for host alignment
HOSTSAM=${HOSTPATH}/${SAMPLE}.${HOSTORG}.bwa.sam
HOSTSAMPE=${HOSTPATH}/${SAMPLE}.${HOSTORG}.bwaPaired.sam


#output files for pAAV2 alignment
VIRALSAM=${VIRALPATH}/${SAMPLE}.${VIRUS}.bwa.sam
VIRALSAMPE=${VIRALPATH}/${SAMPLE}.${VIRUS}.bwaPaired.sam

echo "host index: ${HOSTINDEX}"
echo "viral index: ${VIRALINDEX}"


#align single reads with BWA
python ./alignSingleReadsWithBWA.py --index $HOSTINDEX --reads $READS --output $HOSTSAM --threshold 10 --hflag 200 --bwa $BWAPATH
python ./alignSingleReadsWithBWA.py --index $VIRALINDEX --reads $READS --output $VIRALSAM --threshold 10 --hflag 200 --bwa $BWAPATH


#aligned paired reads with BWA
python ./alignPEReadsWithBWA.py --index $HOSTINDEX --read1 $READ1 --read2 $READ2 --output $HOSTSAMPE --threshold 10 --hflag 200 --bwa $BWAPATH
python ./alignPEReadsWithBWA.py --index $VIRALINDEX --read1 $READ1 --read2 $READ2 --output $VIRALSAMPE --threshold 10 --hflag 200 --bwa $BWAPATH

