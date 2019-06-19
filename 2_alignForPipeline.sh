#!/bin/bash

## Perform alignments to host and viral genomes for pipeline analysis


#SBATCH --job-name=align			#name of job
#SBATCH --time=4:00:00  				#time allowed
#SBATCH --mem-per-cpu=5GB			#memory per cpu
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=3
#SBATCH --mail-type=END,FAIL			#Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=suzanne.scott@csiro.au
#SBATCH --output=../slogs/%x_%j.log 		#Standard error and output log

pwd; hostname; date

#modules
module load parallel

#input and output directories
PROJ=/datastore/sco305/integration/expt2_pipeline-tweaks
DATA=${PROJ}/data
BWAPATH=${PROJ}/tools/bwa-0.7.17/bwa

#get directories with read data
#CSV with data directory in first position, sample name in second, host genome in third and virus genome in fourth
echo -n "" > ${DATA}/data2analyse.txt
dirs=($(find ${DATA}/* -maxdepth 0 -type d))
for dir in "${dirs[@]}"; do
  mkdir -p ${dir}/aligned
  while IFS= read -r line; do
	echo "${dir},${line}" >> ${DATA}/data2analyse.txt
  done < "${dir}/samples.txt"
done

TOTALTASKS=$((${SLURM_NTASKS_PER_NODE} * ${SLURM_NNODES}))
echo total tasks: $TOTALTASKS

# This specifies the options used to run srun. The "-N1 -n${SLURM_NTASKS_PER_NODE} -c1" options are
# used to allocate ${SLURM_NTASKS_PER_NODE} tasks and 1 cpu per task to each node.
srun="srun --exclusive -N1 -n1"

# This specifies the options used to run GNU parallel:
#
#   --delay of 0.2 prevents overloading the controlling node.
#
#   -j is the number of tasks run simultaneously.
#
#   The combination of --joblog and --resume create a task log that
#   can be used to monitor progress.
#
parallel="parallel --delay 0.2 -j $TOTALTASKS --joblog ../slogs/runtask_${SLURM_JOB_ID}.log --resume"

# Run a script, do_merging.sh, using GNU parallel and srun. Parallel
# will run the runtask script for the numbers 1 through ${NAACS}. To
# illustrate, the first job will run like this:
#
#   srun --exclusive -N1 -n1 ./do_merging.sh ${DATA}/firstdataset  > ${DATA}/firstdataset/merge.log
#
$parallel "$srun ./do_aligning.sh {} ${PROJ} ${BWAPATH}  2>&1 {}/align.log" :::: ${DATA}/data2analyse.txt

wait 


date
