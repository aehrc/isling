#!/bin/bash

# This script runs pipeline on alignments run with different parameters
# Pass in path to data as first argument, accession number as second argument, project path as third input 

#modules
module load bedtools
module load samtools



#input paths
DATA=${1}
SAMPLE=$2
SAMPLIST=${1}/samples.txt
OUTPATH=${1}/out/htest
PROJ=$3

#get host and viral references for alignment
HOST=(`grep ${SAMPLE} ${SAMPLIST} | cut -f2 -d,`)
VIRUS=(`grep ${SAMPLE} ${SAMPLIST} | cut -f3 -d,`)

HOSTSAM=${DATA}/hostAligned/${SAMPLE}.${HOST}.bwa.sam
VIRALSAM=${DATA}/viralAligned/${SAMPLE}.${VIRUS}.bwa.sam

HOSTSAM20=${DATA}/hostAligned/${SAMPLE}.${HOST}.bwa.h20.sam
VIRALSAM20=${DATA}/viralAligned/${SAMPLE}.${VIRUS}.bwa.h20.bwa.sam

HOSTSAM10=${DATA}/hostAligned/${SAMPLE}.${HOST}.bwa.h10.sam
VIRALSAM10=${DATA}/viralAligned/${SAMPLE}.${VIRUS}.bwa.h10.sam

echo analysing sample $SAMPLE after alignment to $HOST and ${VIRUS} genomes


#run pipeline
echo ""
echo running integration site anaylsis script on ${SAMPLE}
perl ${PROJ}/tools/AAVintegrationPipeline_master.pl --viral ${VIRALSAM} --human ${HOSTSAM} --output ${OUTPATH}/${SAMPLE}.integrations.txt --bed ${OUTPATH}/${SAMPLE}.integrations.bed --merged ${OUTPATH}/${SAMPLE}.integrations.merged.bed --cutoff 10

perl ${PROJ}/tools/AAVintegrationPipeline_master.pl --viral ${VIRALSAM20} --human ${HOSTSAM20} --output ${OUTPATH}/${SAMPLE}.integrations.h20.txt --bed ${OUTPATH}/${SAMPLE}.integrations.h20.bed --merged ${OUTPATH}/${SAMPLE}.integrations.merged.h1.bed --cutoff 10

perl ${PROJ}/tools/AAVintegrationPipeline_master.pl --viral ${VIRALSAM10} --human ${HOSTSAM10} --output ${OUTPATH}/${SAMPLE}.integrations.h10.txt --bed ${OUTPATH}/${SAMPLE}.integrations.h10.bed --merged ${OUTPATH}/${SAMPLE}.integrations.merged.h10.bed --cutoff 10

sed -i 's\chr\\g' ${OUTPATH}/${SAMPLE}.integrations.txt
sed -i 's\chr\\g' ${OUTPATH}/${SAMPLE}.integrations.h1.txt
sed -i 's\chr\\g' ${OUTPATH}/${SAMPLE}.integrations.h10.txt

#sort output
sort -n -k1,1 -k2,2 ${OUTPATH}/${SAMPLE}.integrations.txt > ${OUTPATH}/${SAMPLE}.integrations.txt.tmp
mv ${OUTPATH}/${SAMPLE}.integrations.txt.tmp ${OUTPATH}/${SAMPLE}.integrations.txt


sort -n -k1,1 -k2,2 ${OUTPATH}/${SAMPLE}.integrations.h1.txt > ${OUTPATH}/${SAMPLE}.integrations.h20.txt.tmp
mv ${OUTPATH}/${SAMPLE}.integrations.h20.txt.tmp ${OUTPATH}/${SAMPLE}.integrations.h20.txt

sort -n -k1,1 -k2,2 ${OUTPATH}/${SAMPLE}.integrations.h10.txt > ${OUTPATH}/${SAMPLE}.integrations.h10.txt.tmp
mv ${OUTPATH}/${SAMPLE}.integrations.h10.txt.tmp ${OUTPATH}/${SAMPLE}.integrations.h10.txt


date
