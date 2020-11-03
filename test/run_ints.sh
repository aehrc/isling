#!/bin/bash

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate snakemake
module load singularity/3.6.4

cd ..

snakemake --configfile test/config/test.yml --cores 1 --rerun-incomplete --use-singularity



#snakemake --configfile ../config/test_pipeline.yml -j 100 --rerun-incomplete --cluster-config cluster.json --use-conda --profile slurm -np

#snakemake --configfile ../config/test_pipeline.yml -j 100 --rerun-incomplete --cluster-config cluster.json --use-conda --profile slurm --dag | dot -Tsvg > dag.svg

#snakemake --configfile ../proj-spcfc/dsets.yaml -j 100 --rerun-incomplete --cluster-config cluster.json --use-conda --profile slurm --dag | dot -Tsvg > dag.svg

#snakemake --configfile ../config/all-cmri.yml -j 100 --rerun-incomplete --cluster-config cluster.json --use-conda --profile slurm -np
