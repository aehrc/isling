#!/bin/bash

# This script runs integration pipeline merged reads to human and viral genomes
# Pass in path to data as first argument, and cpus_per_task parameter from parent sbatch job as second parameter


pwd ; hostname ; date

#load modules
module load parallel

DATA=$1
SLURM_CPUS_PER_TASK=$2
PROJ=$3

#get sample info
SAMPLIST=${DATA}/samples.txt

#make array from columns of SAMPLIST FILE
mapfile SAMPS < <(cut -f1 -d, ${SAMPLIST})

#number of directories
NDIR=(${#SAMPS[@]})

mkdir -p ${DATA}/out



# This specifies the options used to run srun. The "-N1 -n1" options are
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
parallel="parallel --delay 0.2 -j $SLURM_CPUS_PER_TASK --joblog ../slogs/runtask_${SLURM_JOB_ID}_${PARALLEL_SEQ}.log --resume"

# Run a script, do_merging.sh, using GNU parallel and srun. Parallel
# will run the runtask script for the numbers 1 through ${NAACS}. To
# illustrate, the first job will run like this:
#
#   srun --exclusive -N1 -n1 ./do_pipeline.sh ${DATA}/firstdataset   > ${DATA}/firstdataset/pipeline.log
#
$parallel "$srun ./do_pipeline.sh ${DATA} {} ${PROJ} > ${DATA}/pipeline_{}.log" ::: ${SAMPLES[@]}
