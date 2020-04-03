#!/bin/bash

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate snakemake

snakemake --configfile $1 -j 100 --rerun-incomplete --cluster-config cluster.json --use-conda --profile slurm --dry-run "${a[@]:1}"

