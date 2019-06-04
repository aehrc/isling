#!/bin/bash

pwd; hostname; date

#need seqprep via conda
source ~/.bashrc
conda activate bioinfo2

DATA=${1}/reads
OUTPATH=${1}/merged_reads
SLURM_CPUS_PER_TASK=$2

mkdir -p ${OUTPATH}

echo "processing ${DATA}"

#get sample info
SAMPLIST=${1}/samples.txt

NSAMPS=$(cat ${SAMPLIST} | wc -l)

#make array from first column of SAMPLIST FILE
mapfile SAMPLES < <(cut -f1 -d, ${SAMPLIST})

echo ${SAMPLES[@]}

#get paths to reads
READ1="_R1.fastq.gz"
READ2="_R2.fastq.gz"

#generate paths for output
OUT1=".seqPrep_processed.R1.fastq.gz"
OUT2=".seqPrep_processed.R2.fastq.gz"
OUT3=".seqPrep_processed.merged.fastq.gz"
OUT=".seqPrep_processed.fastq"

# This specifies the options used to run srun. The "-N1 -n1 -c1" options are
# used to allocates a single core to each task.
srun="srun --exclusive -N1 -n1 -c1"

# This specifies the options used to run GNU parallel:
#
#   --delay of 0.2 prevents overloading the controlling node.
#
#   -j is the number of tasks run simultaneously.
#
#   The combination of --joblog and --resume create a task log that
#   can be used to monitor progress.
#
parallel1="parallel --delay 0.2 -j $SLURM_CPUS_PER_TASK --joblog ../slogs/runtask_${SLURM_JOB_ID}_${PARALLEL_SEQ}_SeqPrep.log --resume"
parallel2="parallel --delay 0.2 -j $SLURM_CPUS_PER_TASK --joblog ../slogs/runtask_${SLURM_JOB_ID}_${PARALLEL_SEQ}_gzip.log --resume"

# Run a script, do_merging.sh, using GNU parallel and srun. Parallel
# will run the runtask script for the numbers 1 through ${NAACS}. To
# illustrate, the first job will run like this:
#
#   srun --exclusive -N1 -n1 ./do_aligning.sh ${DATA}/firstdataset  > ${OUTPATH}/accession.log
#
$parallel1 "$srun SeqPrep -f ${DATA}/{}${READ1} -r ${DATA}/{}${READ2} -1 ${OUTPATH}/{}${OUT1} -2 ${OUTPATH}/{}${OUT2} -s ${OUTPATH}/{}${OUT3}" ::: ${SAMPLES[@]}

wait

$parallel2 "$srun gzip -cdk ${OUTPATH}/{}* > ${OUTPATH}/{}${OUT}" ::: ${SAMPLES[@]}


date

