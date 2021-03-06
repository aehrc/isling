#!/bin/bash

# dependencies:
## conda
## singularity

set -euo pipefail

CORES="15"

#### dependencies ####

# create conda environment for downloads and snakemake, if it doesn't already exist
conda list -n snakemake_sra || mamba env create -f conda/sra.yml

eval "$(conda shell.bash hook)"
conda activate snakemake_sra

#### references ####

# download large references
echo "downloading references"
bash src/download_references/download_refs.sh

#### reads ####

# download for PRJNA485509 (Nelson)
echo "downloading reads for PRJNA485509"
bash src/download_sra/run_get_SRA.sh config/download/PRJNA485509.yml $CORES

# download for SRP023539" (Lau)
echo "downloading reads for SRP023539"
bash src/download_sra/download_SRP023539.sh $CORES

# download for PRJEB2869 (Sung)
echo "downloading reads for PRJEB2869"
bash src/download_sra/run_get_SRA.sh config/download/PRJEB2869.yml $CORES

# download for PRJNA606282 (Ngyuen)
echo "downloading reads for PRJNA606282"
bash src/download_sra/get_ENA.sh \
data/metadata/PRJNA606282/filereport_read_run_PRJNA606282_tsv.txt \
data/reads/PRJNA606282 \
pearcey-login


#### isling ####

# run PRJNA485509 (Nelson)
echo "running isling on PRJNA485509"
bash src/run_isling/run_isling.sh benchmark/sra_data/config/isling/Nelson_2019_PRJNA485509.yml $CORES

# run SRP023539 (Lau)
echo "running isling on SRP023539"
bash src/run_isling/run_isling.sh benchmark/sra_data/config/isling/Lau_2014_SRP023539_hg19.yml $CORES 

# run PRJEB2869 (Sung)
echo "running isling on PRJEB2869"
bash src/run_isling/run_isling.sh benchmark/sra_data/config/isling/Sung_2012_PRJEB2869.yml $CORES

# run PRJNA606282 (Ngyuen)
echo "running isling on PRJNA606282"
bash src/run_isling/run_isling.sh benchmark/sra_data/config/isling/Ngyuen_2020_PRJNA606282_seed14.yml $CORES

#### other tools ####

echo "running other viral integration tools on PRJNA485509"
bash src/run_other_tools/run_other_tools.sh config/other_tools/Nelson_2019_PRJNA485509.yml $CORES

# run SRP023539 (Lau)
echo "running other viral integration tools on SRP023539"
bash src/run_other_tools/run_other_tools.sh config/other_tools/Lau_2014_SRP023539_hg19.yml $CORES

# run SRP023539 (Lau)
echo "running other viral integration tools on PRJEB2869"
bash src/run_other_tools/run_other_tools.sh config/other_tools/Sung_2012_PRJEB2869.yml $CORES

# run PRJNA606282 (Ngyuen)
echo "running isling on PRJNA606282"
bash src/run_other_tools/run_other_tools.sh config/other_tools/Ngyuen_2020_PRJNA606282.yml $CORES

#### generate figures ####

echo "generating figures and tables"
bash src/analysis/run_Rscript.sh

