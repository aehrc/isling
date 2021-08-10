#!/bin/bash

set -euo pipefail

CONTAINER="rscripts"
VERSION="5"

if [ ! -e ${CONTAINER}_${VERSION}.sif ] ; then
	singularity pull docker://szsctt/${CONTAINER}:${VERSION}
fi

cd src/analysis

singularity exec -B "$(realpath ../..)" ../../${CONTAINER}_${VERSION}.sif  Rscript make_figures_and_tables.R
