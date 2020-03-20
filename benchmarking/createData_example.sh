#!/bin/bash
#
#SBATCH --job-name=simulation		#name of job
#SBATCH --mail-type=BEGIN,END,FAIL		#email the events (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --output=output_location.txt 		#standard error and output log
#SBATCH --mail-user=youremail@csiro.au		#where to send mail
#SBATCH --nodes=1				#number of nodes 
#SBATCH --ntasks=1				#number of tasks 
#SBATCH --cpus-per-task=1			#number of CPUs per task
#SBATCH	--mem=12GB 				#memory per NODE 
#SBATCH	--time=01:30:00				#time limit hrs:min:sec
#SBATCH	--error=error_location.txt 		#error destination 

pwd; hostname; date

echo "Running on $SLURM_JOB_NUM_NODES nodes with $SLURM_NTASKS tasks, each with $SLURM_CPUS_PER_TASK cores."

echo "Batch is $SBATCH_BATCH, job is $SBATCH_JOB_NAME"


srun -N 1 python3 insert_virus.py --host test_host.fa --virus test_virus.fa --ints test_ints.fa --locs test_locs.txt --host_ints test_host_ints.csv --int_num 1 

srun -N 1 art_illumina -ss HS25 -sam -i test_ints.fa -p -l 150 -f 20 -m 400 -s 15 -o test_data

srun -N 1 python3 read_integrations.py --sam test_data.sam --host_ints test_host_ints.csv --save all_reads.csv --viral_reads viral_reads.csv

wait

date
