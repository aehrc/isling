#!/bin/bash
set -euo pipefail

# usage: run_isling.sh <CONFIGFILE>

CONFIG="$1"

cd ../..

snakemake --keep-going --rerun-incomplete --use-singularity --jobs 1 --configfile $CONFIG
