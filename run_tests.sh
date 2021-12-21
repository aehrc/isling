#!/bin/bash

# run test for Docker container
snakemake \
	-j 1 \
	--forceall \
	--configfile test/config/test.yml

# container has micromamba, not mamba or conda, so we need a workaround
export PATH="/opt/conda/bin/conda:$PATH" 

snakemake \
	-j 1 \
	--forceall \
	--configfile test/config/test.yml \
	--use-conda \
	--conda-frontend conda
