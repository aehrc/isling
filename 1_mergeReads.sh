#!/bin/bash

## process and merge overlapping reads with SeqPrep
## Compress to a single fasta file

#SBATCH --job-name=merge-for-LWpipeline		#name of job
#SBATCH --time=12:00:00  			#time allowed
#SBATCH --mem=40GB				#memory
#SBATCH --nodes=1				#number of nodes
#SBATCH --ntasks-per-node=1			#number of tasks per node
#SBATCH -a 1-12					#sbatch array - number of tasks
#SBATCH --mail-type=END,FAIL			#Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=suzanne.scott@csiro.au
#SBATCH --output=../slogs/%x_%j.log 		#Standard error and output log

pwd; hostname; date

#need seqprep via conda
source ~/.bashrc
conda activate bioinfo2

#input and output directories
PROJ=/datastore/sco305/integration/expt2_testing-secondary-alignments
DATA=/datastore/sco305/integration/expt1_comparing-pipelines/data/tool_comparison/rawdata
OUTPATH=${PROJ}/preproc/merged

mkdir -p $OUTPATH
cp ${DATA}/sample2index.txt ${PROJ}/preproc
SUFFIX=$(head /dev/urandom | tr -dc a-z | head -c3)
sed 's/human/hg38/g; s/hg19/hg38/g; s/references/expt2_testing-secondary-alignments\/references/g' ${PROJ}/preproc/sample2index.txt > ${PROJ}/preproc/sample2index.txt.${SUFFIX}

#get sample names
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${PROJ}/preproc/sample2index.txt.${SUFFIX} | cut -d',' -f1)


echo $SAMPLE

#get paths to reads
READ1=${DATA}/${SAMPLE}_R1.fastq.gz
READ2=${DATA}/${SAMPLE}_R2.fastq.gz

#generate paths for output
OUT1=${OUTPATH}/${SAMPLE}.seqPrep_processed.R1.fastq.gz
OUT2=${OUTPATH}/${SAMPLE}.seqPrep_processed.R2.fastq.gz
OUT3=${OUTPATH}/${SAMPLE}.seqPrep_processed.merged.fastq.gz

OUT=${OUTPATH}/${SAMPLE}.seqPrep_processed.fastq

#run seqprep
SeqPrep -f $READ1 -r $READ2 -1 $OUT1 -2 $OUT2 -s $OUT3

gzip -cdk ${OUTPATH}/${SAMPLE}* > $OUT

mv ${PROJ}/preproc/sample2index.txt.${SUFFIX} ${PROJ}/preproc/sample2index.txt

date
