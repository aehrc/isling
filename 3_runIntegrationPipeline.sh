#!/bin/bash

## Runs the AAV integration pipeline on alignment files

#SBATCH --job-name=run-LWpipeline		#name of job
#SBATCH --time=0:30:00  			#time allowed
#SBATCH --mem=10GB				#memory
#SBATCH --nodes=1				#number of nodes
#SBATCH --ntasks-per-node=1			#number of tasks per node
#SBATCH -a 1-12					#sbatch array - number of tasks
#SBATCH --mail-type=END,FAIL			#Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=suzanne.scott@csiro.au
#SBATCH --output=../slogs/%x_%j.log 		#Standard error and output log

pwd; hostname; date

module load bedtools
module load samtools

PROJ=/datastore/sco305/integration/expt2_testing-secondary-alignments
DATA=${PROJ}/preproc

SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${DATA}/sample2index.txt | cut -d',' -f1)
HOST=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${DATA}/sample2index.txt | cut -d',' -f3)

echo $SAMPLE

HOSTSAM=${DATA}/hostAligned/${SAMPLE}.${HOST}.bwa.sorted.supFilt.sam
VIRALSAM=${DATA}/viralAligned/${SAMPLE}.virus.bwa.sorted.supFilt.sam
AHOSTSAM=${DATA}/hostAligned/${SAMPLE}.${HOST}.bwa.a.sorted.supFilt.sam
AVIRALSAM=${DATA}/viralAligned/${SAMPLE}.virus.bwa.a.sorted.supFilt.sam
HHOSTSAM=${DATA}/hostAligned/${SAMPLE}.${HOST}.bwa.h.sorted.supFilt.sam
HVIRALSAM=${DATA}/viralAligned/${SAMPLE}.virus.bwa.h.sorted.supFilt.sam

OUTPATH=${PROJ}/out

mkdir -p $OUTPATH

#run pipeline
echo ""
echo running normal alignment
perl ${PROJ}/tools/AAVintegrationPipeline_master.pl --viral ${VIRALSAM} --human ${HOSTSAM} --output ${OUTPATH}/${SAMPLE}.integrations.txt --bed ${OUTPATH}/${SAMPLE}.integrations.bed --merged ${OUTPATH}/${SAMPLE}.integrations.merged.bed --cutoff 10
echo ""
echo running with a flag
perl ${PROJ}/tools/AAVintegrationPipeline_master.pl --viral ${AVIRALSAM} --human ${AHOSTSAM} --output ${OUTPATH}/${SAMPLE}.integrations.a.txt --bed ${OUTPATH}/${SAMPLE}.integrations.a.bed --merged ${OUTPATH}/${SAMPLE}.integrations.a.merged.bed --cutoff 10
echo ""
echo running with h flag
perl ${PROJ}/tools/AAVintegrationPipeline_master.pl --viral ${HVIRALSAM} --human ${HHOSTSAM} --output ${OUTPATH}/${SAMPLE}.integrations.h.txt --bed ${OUTPATH}/${SAMPLE}.integrations.h.bed --merged ${OUTPATH}/${SAMPLE}.integrations.merged.h.bed --cutoff 10

#sort output
sort -k1,1 -k2,2 ${OUTPATH}/${SAMPLE}.integrations.txt > ${OUTPATH}/${SAMPLE}.integrations.txt
sort -k1,1 -k2,2 ${OUTPATH}/${SAMPLE}.integrations.a.txt > ${OUTPATH}/${SAMPLE}.integrations.a.txt
sort -k1,1 -k2,2 ${OUTPATH}/${SAMPLE}.integrations.h.txt > ${OUTPATH}/${SAMPLE}.integrations.h.txt

#get read IDs
cat ${OUTPATH}/${SAMPLE}.integrations.txt | cut -f12 | tail -n+1 > ${OUTPATH}/${SAMPLE}.readIDs.txt

echo extracting reads...

#get fastq of reads
samtools view -H ${HOSTSAM} > ${OUTPATH}/${SAMPLE}.reads.sam
samtools view ${HOSTSAM} | fgrep -w -f ${OUTPATH}/${SAMPLE}.readIDs.txt >> ${OUTPATH}/${SAMPLE}.reads.sam

samtools fastq ${OUTPATH}/${SAMPLE}.reads.sam > ${OUTPATH}/${SAMPLE}.readsofInterest.fastq

#rm ${OUTPATH}/${SAMPLE}.reads.sam ${OUTPATH}/${SAMPLE}.readIDs.txt

date
