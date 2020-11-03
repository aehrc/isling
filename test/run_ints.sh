#!/bin/bash

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate snakemake

snakemake --configfile $1 -j 100 --rerun-incomplete --cluster-config cluster.json --use-conda --profile slurm "${a[@]:1}"

#snakemake --configfile ../config/test_pipeline.yml -j 100 --rerun-incomplete --cluster-config cluster.json --use-conda --profile slurm -np

#snakemake --configfile ../config/test_pipeline.yml -j 100 --rerun-incomplete --cluster-config cluster.json --use-conda --profile slurm --dag | dot -Tsvg > dag.svg

#snakemake --configfile ../proj-spcfc/dsets.yaml -j 100 --rerun-incomplete --cluster-config cluster.json --use-conda --profile slurm --dag | dot -Tsvg > dag.svg

#snakemake --configfile ../config/all-cmri.yml -j 100 --rerun-incomplete --cluster-config cluster.json --use-conda --profile slurm -np
