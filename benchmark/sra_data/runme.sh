#!/bin/bash

# dependencies:
## conda
## singularity
## parallel

set -euo pipefail

module load parallel
module load singularity

#### references ####

# download large references
echo "downloading references"
bash src/download_references/download_refs.sh

#### reads ####

# create conda environment for downloads, if it doesn't already exist
conda list -n snakemake_sra || conda env create -f conda/sra.yml

eval "$(conda shell.bash hook)"
conda activate snakemake_sra

# download for PRJNA485509 (Nelson)
echo "downloading reads for PRJNA485509"
bash src/download_sra/run_get_SRA.sh config/download/PRJNA485509.yml

# download for SRP023539" (Lau)
echo "downloading reads for SRP023539"
bash src/download_sra/download_SRP023539.sh

# download for PRJEB2869 (Sung)
echo "downloading reads for PRJEB2869"
bash src/download_sra/run_get_SRA.sh config/download/PRJEB2869.yml

#### isling ####

# run PRJNA485509 (Nelson)
echo "running isling on PRJNA485509"
bash src/run_isling/run_isling.sh benchmark/sra_data/config/isling/Nelson_2019_PRJNA485509.yml

