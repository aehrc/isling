#!/bin/bash

#modules
module load bedtools
module load samtools

#get info from inputs
DATADIR=$(echo $1 | cut -f1 -d,)
DATA="${DATADIR}/aligned"
OUTPATH="${DATADIR}/out"
SAMPLE=$(echo $1 | cut -f2 -d,)
HOST=$(echo $1 | cut -f3 -d,)
VIRUS=$(echo $1 | cut -f4 -d,)
PROJ=$2

HOSTSAM=${DATA}/${SAMPLE}.${HOST}.bwa.sam
VIRALSAM=${DATA}/${SAMPLE}.${VIRUS}.bwa.sam

echo analysing sample $SAMPLE after alignment to ${HOST} and ${VIRUS} genomes


#run pipeline
echo ""
echo running integration site anaylsis script on ${SAMPLE}
perl ${PROJ}/tools/AAVintegrationPipeline_master.pl --viral ${VIRALSAM} --human ${HOSTSAM} --output ${OUTPATH}/${SAMPLE}.${HOST}_${VIRUS}.integrations.txt --bed ${OUTPATH}/${SAMPLE}.${HOST}_${VIRUS}.integrations.bed --merged ${OUTPATH}/${SAMPLE}.${HOST}_${VIRUS}.integrations.merged.bed --cutoff 10

sed -i 's\chr\\g' ${OUTPATH}/${SAMPLE}.${HOST}_${VIRUS}.integrations.txt

sort output
sort -n -k1,1 -k2,2 ${OUTPATH}/${SAMPLE}.${HOST}_${VIRUS}.integrations.txt > ${OUTPATH}/${SAMPLE}.${HOST}_${VIRUS}.integrations.txt.tmp
mv ${OUTPATH}/${SAMPLE}.${HOST}_${VIRUS}.integrations.txt.tmp ${OUTPATH}/${SAMPLE}.${HOST}_${VIRUS}.integrations.txt



date
