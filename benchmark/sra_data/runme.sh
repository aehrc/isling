#!/bin/bash

set -euo pipefail

# download large references
echo "downloading references"
bash src/download_references/download_refs.sh

# create conda environment for downloads, if it doesn't already exist
conda list -n snakemake_sra || conda env create -f conda/sra.yml

# download for PRJNA485509 (Nelson)
echo "downloading reads for PRJNA485509"
bash src/download_sra/run_get_SRA.sh config/download/PRJNA485509.yml
