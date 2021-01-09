#!/bin/bash

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate snakemake

module load singularity

cd ..

snakemake \
	-j 1 \
	--rerun-incomplete \
	--use-singularity \
	--configfile test/config/test.yml
