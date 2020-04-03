#!/bin/bash

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate snakemake

snakemake --configfile $1 -j 100 -s viralIntegrations.sf --rerun-incomplete --cluster-config cluster.json --use-conda --profile slurm "${a[@]:1}"

# --restart-times 3
