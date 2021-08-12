#!/bin/bash
set -euo pipefail

# usage: run_isling.sh <CONFIGFILE> <CORES>

CONFIG="$1"
CORES="$2"

cd ../..

snakemake --keep-going --rerun-incomplete --use-singularity --cores $CORES --configfile $CONFIG
