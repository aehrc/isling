#!/bin/bash

# This script aligns merged reads to human and viral genomes, testing different paraemeters for alignment
# Pass in path to data as first argument, accession number as second argument, project path as third input 

pwd ; hostname ; date

#load modules
module load samtools
module load bwa

#need conda for pysam
source ~/.bashrc
conda activate bioinfo2

#input paths
DATA=${1}/merged_reads
SAMPLE=$2
SAMPLIST=${1}/samples.txt
PROJ=${3}

#location of reads - here only using merged reads
READS=${DATA}/${SAMPLE}.seqPrep_processed.merged.fastq.gz

#output paths
HOSTPATH=${1}/hostAligned
VIRALPATH=${1}/viralAligned

echo data path: ${DATA}
echo sample: ${SAMPLE}
echo sample list: ${SAMPLIST}
echo project directory: ${PROJ}


#get host and viral references for alignment
HOST=(`grep ${SAMPLE} ${SAMPLIST} | cut -f2 -d,`)
VIRUS=(`grep ${SAMPLE} ${SAMPLIST} | cut -f3 -d,`)

#index locations
HOSTINDEX=${PROJ}/references/${HOST}
VIRALINDEX=${PROJ}/references/${VIRUS}


#output files for host alignment (regular alignment)
HOSTSAM=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.sam
HOSTSAM1=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.h20.sam
HOSTSAM10=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.h10.sam


#output files for pAAV2 alignment
VIRALSAM=${VIRALPATH}/${SAMPLE}.${VIRUS}.bwa.sam
VIRALSAM1=${VIRALPATH}/${SAMPLE}.${VIRUS}.h20.bwa.sam
VIRALSAM10=${VIRALPATH}/${SAMPLE}.${VIRUS}.h10.bwa.sam


#align single reads with BWA
#regular alignment
echo aligning...
python ./alignSingleReadsWithBWA.py --index $HOSTINDEX --reads $READS --output $HOSTSAM --threshold 10 
python ./alignSingleReadsWithBWA.py --index $VIRALINDEX --reads $READS --output $VIRALSAM --threshold 10 

python ./alignSingleReadsWithBWA.py --index $HOSTINDEX --reads $READS --output $HOSTSAM1 --threshold 10 --hflag 20
python ./alignSingleReadsWithBWA.py --index $VIRALINDEX --reads $READS --output $VIRALSAM1 --threshold 10 --hflag 20

python ./alignSingleReadsWithBWA.py --index $HOSTINDEX --reads $READS --output $HOSTSAM10 --threshold 10 --hflag 10
python ./alignSingleReadsWithBWA.py --index $VIRALINDEX --reads $READS --output $VIRALSAM10 --threshold 10 --hflag 10

