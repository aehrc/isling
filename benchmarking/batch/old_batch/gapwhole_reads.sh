#!/bin/bash
#
#SBATCH --job-name=integrate_test		#name of job
#SBATCH --mail-type=BEGIN,END,FAIL		#email the events (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --output=../../slogs/chr20/gapwhole/300int_rep3.txt #standard error and output log
#SBATCH --mail-user=Susie.Grigson@cisro.au	#where to send mail
#SBATCH --nodes=1				#number of nodes 
#SBATCH --ntasks=1				#number of tasks 
#SBATCH --cpus-per-task=1			#number of CPUs per task
#SBATCH	--mem=12GB 				#memory per NODE 
#SBATCH	--time=00:30:00				#time limit hrs:min:sec
#SBATCH	--error=../../slogs/chr20/gapwhole/300int_rep3_error.txt 		#error destination 

pwd; hostname; date

echo "Running on $SLURM_JOB_NUM_NODES nodes with $SLURM_NTASKS tasks, each with $SLURM_CPUS_PER_TASK cores."

echo "Batch is $SBATCH_BATCH, job is $SBATCH_JOB_NAME"

srun -N 1 python3 insert_virus.py --host lower_chr20.fna --virus AAV_vector.fa --ints ../../script_outputs/chr20/gapwhole/300int/rep_3/ints.fa --locs ../../script_outputs/chr20/gapwhole/300int/rep_3/locs.txt --ints_host ../../script_outputs/chr20/gapwhole/300int/rep_3/host_integrations.csv --int_num 300 
srun -N 1 art_illumina -ss HS25 -sam -i ../../script_outputs/chr20/gapwhole/300int/rep_3/ints.fa -p -l 150 -f 20 -m 400 -s 15 -o ../../script_outputs/chr20/gapwhole/300int/rep_3/paired_dat
