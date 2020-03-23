#!/bin/bash

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate snakemake

snakemake -j 1 -s viralIntegrations.sf --rerun-incomplete --use-conda --configfile ../proj-spcfc/test.yml
