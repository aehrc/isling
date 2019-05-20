#!/bin/bash

## index references

#SBATCH --job-name=index-refs			#name of job
#SBATCH --time=12:00:00  			#time allowed
#SBATCH --mem=40GB				#memory
#SBATCH --nodes=3				#number of nodes
#SBATCh --ntasks=3
#SBATCH --ntasks-per-node=1			#number of tasks per node
#SBATCH --mail-type=END,FAIL			#Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=suzanne.scott@csiro.au
#SBATCH --output=../slogs/%x_%j.log 		#Standard error and output log

pwd; hostname; date

module load bwa

#input and output directories
PROJ=/datastore/sco305/integration/expt2_testing-secondary-alignments
REFS=${PROJ}/references

srun -N 1 -n 1 bwa index ${REFS}/mouse.fa &
srun -N 1 -n 1 bwa index ${REFS}/hg38.fa &
srun -N 1 -n 1 bwa index ${REFS}/virus.fa

wait

date
