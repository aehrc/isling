#!/bin/bash


# This script aligns merged reads to human and viral genomes
# Pass in (path to data, sample name, host reference, virus reference) as first argument, project directory as second argument

pwd ; hostname ; date

#get info from inputs
DATADIR=$(echo $1 | cut -f1 -d,)
DATA="${DATADIR}/merged_reads"
SAMPLE=$(echo $1 | cut -f2 -d,)
HOST=$(echo $1 | cut -f3 -d,)
VIRUS=$(echo $1 | cut -f4 -d,)
PROJ=$2

#load modules
module load samtools
module load bwa

#need conda for pysam
source ~/.bashrc
conda activate bioinfo2

#input paths


#location of reads - here only using merged reads
READS=${DATA}/${SAMPLE}.seqPrep_processed.fastq

#output paths
HOSTPATH=${DATADIR}/aligned
VIRALPATH=${DATADIR}/aligned


#index locations
HOSTINDEX=${PROJ}/references/${HOST}
VIRALINDEX=${PROJ}/references/${VIRUS}


#output files for host alignment (regular alignment)
HOSTSAM=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.sam
HOSTBAM=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.bam
HOSTSOR=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.sorted.bam
HOSTSUP=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.sorted.supFilt.bam
HOSTSUPSAM=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.sorted.supFilt.sam


#output files for pAAV2 alignment
VIRALSAM=${VIRALPATH}/${SAMPLE}.${VIRUS}.bwa.sam
VIRALBAM=${VIRALPATH}/${SAMPLE}.${VIRUS}.bwa.bam
VIRALSOR=${VIRALPATH}/${SAMPLE}.${VIRUS}.bwa.sorted.bam
VIRALSUP=${VIRALPATH}/${SAMPLE}.${VIRUS}.bwa.sorted.supFilt.bam
VIRALSUPSAM=${VIRALPATH}/${SAMPLE}.${VIRUS}.bwa.sorted.supFilt.sam


#align single reads with BWA
python ./alignSingleReadsWithBWA.py --index $HOSTINDEX --reads $READS --output $HOSTSAM --threshold 10 --hflag 20 
python ./alignSingleReadsWithBWA.py --index $VIRALINDEX --reads $READS --output $VIRALSAM --threshold 10 --hflag 20


#convert to bam and sort
python ./convertAndSort.py --sam $HOSTSAM --bam $HOSTBAM --sort $HOSTSOR
python ./convertAndSort.py --sam $VIRALSAM --bam $VIRALBAM --sort $VIRALSOR


#remove supplementary alignments
samtools view -bh -G 0x800 $HOSTSOR > $HOSTSUP
samtools view -bh -G 0x800 $VIRALSOR > $VIRALSUP


#index 
samtools index $HOSTSUP
samtools index $VIRALSUP


#convert to sam
samtools view -h -o $HOSTSUPSAM $HOSTSUP
samtools view -h -o $VIRALSUPSAM $VIRALSUP
