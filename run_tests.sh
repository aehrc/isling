#!/bin/bash

# run test for Docker container
snakemake \
	-j 1 \
	--forceall \
	--configfile test/config/test.yml

snakemake \
	-j 1 \
	--forceall \
	--configfile test/config/test.yml \
	--use-conda
