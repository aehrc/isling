#!/bin/bash

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate snakemake

snakemake --configfile $1 -j 100 --restart-times 3 -s viralIntegrations.sf --rerun-incomplete --use-conda --cluster "sbatch --mem-per-cpu {cluster.mem-per-cpu}  -t {cluster.time} --mail-user {cluster.mail-user} --nodes {cluster.nodes} --cpus-per-task {cluster.cpus-per-task} --output {cluster.output}" $2
