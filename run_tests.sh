#!/bin/bash

# run test for Docker container
snakemake \
	-j 1 \
	--forceall \
	--configfile test/config/test.yml

# container has micromamba, not mamba or conda, so we need a workaround
shopt -s expand_aliases
alias conda=/opt/conda/bin/conda

snakemake \
	-j 1 \
	--forceall \
	--configfile test/config/test.yml \
	--use-conda \
	--conda-frontend conda
