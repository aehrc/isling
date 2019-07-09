#!/bin/bash

#modules
module load bedtools
module load samtools

#get info from inputs
DATADIR=$(echo $1 | cut -f1 -d,)
DATA="${DATADIR}/aligned"
OUTPATH="${DATADIR}/out"
SAMPLE=$(echo $1 | cut -f2 -d,)
HOSTORG=$(echo $1 | cut -f3 -d,)
VIRUS=$(echo $1 | cut -f4 -d,)
PROJ=$2

HOSTSAM=${DATA}/${SAMPLE}.${HOSTORG}.bwa.sam
VIRALSAM=${DATA}/${SAMPLE}.${VIRUS}.bwa.sam
HOSTPESAM=${DATA}/${SAMPLE}.${HOSTORG}.bwaPaired.sam
VIRALPESAM=${DATA}/${SAMPLE}.${VIRUS}.bwaPaired.sam

echo analysing sample $SAMPLE after alignment to ${HOSTORG} and ${VIRUS} genomes

#run pipeline
echo ""
echo running integration site anaylsis script on ${SAMPLE}
perl -I ${PROJ}/tools  ${PROJ}/tools/softClip.pl --viral ${VIRALSAM} --human ${HOSTSAM} --output ${OUTPATH}/${SAMPLE}.${HOSTORG}_${VIRUS}.soft.txt --cutoff 10

perl -I ${PROJ}/tools ${PROJ}/tools/discordant.pl --viral ${VIRALPESAM} --human ${HOSTPESAM} --output ${OUTPATH}/${SAMPLE}.${HOSTORG}_${VIRUS}.discordant.txt 

 

#concatenate output
awk 'FNR==1 && NR!=1 { getline; } { print $0; }' ${OUTPATH}/${SAMPLE}.${HOSTORG}_${VIRUS}*.txt > ${OUTPATH}/${SAMPLE}.${HOSTORG}_${VIRUS}.integrations.txt

sed -i 's\chr\\g' ${OUTPATH}/${SAMPLE}.${HOSTORG}_${VIRUS}.integrations.txt

#sort output
sort -n -k1,1 -k2,2 ${OUTPATH}/${SAMPLE}.${HOSTORG}_${VIRUS}.integrations.txt > ${OUTPATH}/${SAMPLE}.${HOSTORG}_${VIRUS}.integrations.txt.tmp
mv ${OUTPATH}/${SAMPLE}.${HOSTORG}_${VIRUS}.integrations.txt.tmp ${OUTPATH}/${SAMPLE}.${HOSTORG}_${VIRUS}.integrations.txt

#check for vector rearrangements
perl -I ${PROJ}/tools  ${PROJ}/tools/checkVecRearrange.pl --viral ${VIRALSAM} --human ${HOSTSAM} --output ${OUTPATH}/${SAMPLE}.${HOSTORG}_${VIRUS}.rearrange.txt --bed ${OUTPATH}/${SAMPLE}.${HOSTORG}_${VIRUS}.rearrange.bed

date
