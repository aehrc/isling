#!/bin/bash

## Perform alignments to host and viral genomes for pipeline analysis

#SBATCH --job-name=align-for-LWpipeline		#name of job
#SBATCH --time=12:00:00  			#time allowed
#SBATCH --mem=40GB				#memory
#SBATCH --nodes=1				#number of nodes
#SBATCH --ntasks-per-node=1			#number of tasks per node
#SBATCH -a 1-12					#sbatch array - number of tasks
#SBATCH --mail-type=END,FAIL			#Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=suzanne.scott@csiro.au
#SBATCH --output=../slogs/%x_%j.log 		#Standard error and output log

pwd; hostname; date

#load modules
module load samtools
module load bwa

#need conda for pysam
source ~/.bashrc
conda activate bioinfo2

#input and output directories
PROJ=/datastore/sco305/integration/expt2_testing-secondary-alignments
DATA=${PROJ}/preproc/merged
OUT=${PROJ}/preproc/

#create output directory
mkdir -p $OUT


#get sample info from sample2index
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${OUT}/sample2index.txt | cut -d',' -f1)
#HOSTINDEX=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${DATA}/sample2index.txt | cut -d',' -f2)
HOST=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${OUT}/sample2index.txt | cut -d',' -f3)

echo $SAMPLE

#indexed genomes
VIRALINDEX=${PROJ}/references/virus
HOSTINDEX=${PROJ}/references/${HOST}

#location of reads
READS=${DATA}/${SAMPLE}.seqPrep_processed.merged.fastq.gz

#output paths
HOSTPATH=${OUT}/hostAligned
VIRALPATH=${OUT}/viralAligned

mkdir -p $HOSTPATH
mkdir -p $VIRALPATH

#output files for host alignment (regular alignment)
HOSTSAM=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.sam
HOSTBAM=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.bam
HOSTSOR=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.sorted.bam
HOSTSUP=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.sorted.supFilt.bam
HOSTSUPSAM=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.sorted.supFilt.sam

#with -a flag (to output secondary alignments)
AHOSTSAM=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.a.sam
AHOSTBAM=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.a.bam
AHOSTSOR=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.a.sorted.bam
AHOSTSUP=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.a.sorted.supFilt.bam
AHOSTSUPSAM=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.a.sorted.supFilt.sam

#with -h flag (to put more secondary alignments in XA
HHOSTSAM=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.h.sam
HHOSTBAM=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.h.bam
HHOSTSOR=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.h.sorted.bam
HHOSTSUP=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.h.sorted.supFilt.bam
HHOSTSUPSAM=${HOSTPATH}/${SAMPLE}.${HOST}.bwa.h.sorted.supFilt.sam


#output files for viral alignment
VIRALSAM=${VIRALPATH}/${SAMPLE}.virus.bwa.sam
VIRALBAM=${VIRALPATH}/${SAMPLE}.virus.bwa.bam
VIRALSOR=${VIRALPATH}/${SAMPLE}.virus.bwa.sorted.bam
VIRALSUP=${VIRALPATH}/${SAMPLE}.virus.bwa.sorted.supFilt.bam
VIRALSUPSAM=${VIRALPATH}/${SAMPLE}.virus.bwa.sorted.supFilt.sam

AVIRALSAM=${VIRALPATH}/${SAMPLE}.virus.bwa.a.sam
AVIRALBAM=${VIRALPATH}/${SAMPLE}.virus.bwa.a.bam
AVIRALSOR=${VIRALPATH}/${SAMPLE}.virus.bwa.a.sorted.bam
AVIRALSUP=${VIRALPATH}/${SAMPLE}.virus.bwa.a.sorted.supFilt.bam
AVIRALSUPSAM=${VIRALPATH}/${SAMPLE}.virus.bwa.a.sorted.supFilt.sam

HVIRALSAM=${VIRALPATH}/${SAMPLE}.virus.bwa.h.sam
HVIRALBAM=${VIRALPATH}/${SAMPLE}.virus.bwa.h.bam
HVIRALSOR=${VIRALPATH}/${SAMPLE}.virus.bwa.h.sorted.bam
HVIRALSUP=${VIRALPATH}/${SAMPLE}.virus.bwa.h.sorted.supFilt.bam
HVIRALSUPSAM=${VIRALPATH}/${SAMPLE}.virus.bwa.h.sorted.supFilt.sam

#align single reads with BWA
#regular alignment
echo regular alignment
python ./alignSingleReadsWithBWA.py --index $HOSTINDEX --reads $READS --output $HOSTSAM --threshold 10 
python ./alignSingleReadsWithBWA.py --index $VIRALINDEX --reads $READS --output $VIRALSAM --threshold 10 

#align with -a flag
echo a flag
python ./alignSingleReadsWithBWA.py --index $HOSTINDEX --reads $READS --output $AHOSTSAM --threshold 10 --aflag True
python ./alignSingleReadsWithBWA.py --index $VIRALINDEX --reads $READS --output $AVIRALSAM --threshold 10 --aflag True

#align with -h flag = 20
echo h flag
python ./alignSingleReadsWithBWA.py --index $HOSTINDEX --reads $READS --output $HHOSTSAM --threshold 10 --hflag 20
python ./alignSingleReadsWithBWA.py --index $VIRALINDEX --reads $READS --output $HVIRALSAM --threshold 10 --hflag 20



#convert to bam and sort
python ./convertAndSort.py --sam $HOSTSAM --bam $HOSTBAM --sort $HOSTSOR
python ./convertAndSort.py --sam $VIRALSAM --bam $VIRALBAM --sort $VIRALSOR
python ./convertAndSort.py --sam $AHOSTSAM --bam $AHOSTBAM --sort $AHOSTSOR
python ./convertAndSort.py --sam $AVIRALSAM --bam $AVIRALBAM --sort $AVIRALSOR
python ./convertAndSort.py --sam $HHOSTSAM --bam $HHOSTBAM --sort $HHOSTSOR
python ./convertAndSort.py --sam $HVIRALSAM --bam $HVIRALBAM --sort $HVIRALSOR

#remove supplementary reads
samtools view -bh -G 0x800 $HOSTSOR > $HOSTSUP
samtools view -bh -G 0x800 $VIRALSOR > $VIRALSUP
samtools view -bh -G 0x800 $AHOSTSOR > $AHOSTSUP
samtools view -bh -G 0x800 $AVIRALSOR > $AVIRALSUP
samtools view -bh -G 0x800 $HHOSTSOR > $HHOSTSUP
samtools view -bh -G 0x800 $HVIRALSOR > $HVIRALSUP

#index 
samtools index $HOSTSUP
samtools index $VIRALSUP
samtools index $AHOSTSUP
samtools index $AVIRALSUP
samtools index $HHOSTSUP
samtools index $HVIRALSUP

#convert to sam
samtools view -h -o $HOSTSUPSAM $HOSTSUP
samtools view -h -o $VIRALSUPSAM $VIRALSUP
samtools view -h -o $AHOSTSUPSAM $AHOSTSUP
samtools view -h -o $AVIRALSUPSAM $AVIRALSUP
samtools view -h -o $HHOSTSUPSAM $HHOSTSUP
samtools view -h -o $HVIRALSUPSAM $HVIRALSUP

date
