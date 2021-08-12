#!/bin/bash

# usage: run_get_SRA.sh <configfile>

set -e

CONFIG=$1
CORES=$2
CONDA="snakemake_sra"

eval "$(conda shell.bash hook)"
conda activate $CONDA

snakemake --snakefile src/download_sra/get_SRA.sf --cores $CORES --configfile $CONFIG -T 3

