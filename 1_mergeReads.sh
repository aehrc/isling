#!/bin/bash

## process and merge overlapping reads with SeqPrep
## Compress to a single fasta file

#SBATCH --job-name=merge			#name of job
#SBATCH --time=2:00:00  				#time allowed
#SBATCH --mem-per-cpu=5GB			#memory
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=5
#SBATCH --mail-type=END,FAIL			#Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=suzanne.scott@csiro.au
#SBATCH --output=../slogs/%x_%j.log 		#Standard error and output log

pwd; hostname; date

#modules
module load parallel

#input and output directories
PROJ=/datastore/sco305/integration/expt2_pipeline-tweaks
DATA=${PROJ}/data

#get directories with read data
echo "getting directories to merge..."
data2analyse=()
while IFS= read -r -d $'\0'; do
	data2analyse+=("${REPLY:0:-6}")
done < <(find ${DATA}/*/reads -maxdepth 0 -type d -print0)

echo ${data2analyse[@]}

#number of directories
NDIR=(${#data2analyse[@]})
echo "Analysing ${NDIR} directories"

# This specifies the options used to run srun. The "-N1 -n${SLURM_NTASKS_PER_NODE} -c1" options are
# used to allocate ${SLURM_NTASKS_PER_NODE} tasks and 1 cpu per task to each node.
srun="srun --exclusive -N1 -n${SLURM_NTASKS_PER_NODE} -c1"

# This specifies the options used to run GNU parallel:
#
#   --delay of 0.2 prevents overloading the controlling node.
#
#   -j is the number of tasks run simultaneously.
#
#   The combination of --joblog and --resume create a task log that
#   can be used to monitor progress.
#
parallel="parallel --delay 0.2 -j ${SLURM_NTASKS_PER_NODE} --joblog ../slogs/runtask_${SLURM_JOB_ID}.log --resume"

# Run a script, do_merging.sh, using GNU parallel and srun. Parallel
# will run the runtask script for the numbers 1 through ${NAACS}. To
# illustrate, the first job will run like this:
#
#   srun --exclusive -N1 -n1 ./do_merging.sh ${DATA}/firstdataset   > ${DATA}/firstdataset/merge.log
#
$parallel "$srun ./do_merging.sh {} $SLURM_CPUS_PER_TASK 2>&1 {}/merge.log" ::: ${data2analyse[@]}


date
