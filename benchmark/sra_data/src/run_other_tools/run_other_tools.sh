#!/bin/bash
set -euo pipefail

# usage: run_other_tools.sh <CONFIGFILE>

CONFIG="$1"

cd src/run_other_tools/

snakemake --keep-going --rerun-incomplete --use-singularity --singularity-args "-B $(realpath ..)" --jobs 1 --configfile ../../$CONFIG

