#!/bin/bash
set -euo pipefail

# usage: run_other_tools.sh <CONFIGFILE> <CORES>

CONFIG="$1"
CORES="$2"

cd src/run_other_tools/

snakemake --keep-going --rerun-incomplete --use-singularity --singularity-args "-B $(realpath ..)" --cores $CORES --configfile ../../$CONFIG

