#!/bin/bash
set -euo pipefail

# $1 - config file
# $2 - name of bucket
# other args - passed to snakemake (no need for --kubernetes or bucket details)

# script to run isling using kubernetes.  Pass through any other arguments passed the this script to snakemake


CONFIG=$1
BUCKET=$2

if [ $# -gt 2 ]; then
	ARGS="${@:3}"
else
	ARGS=""
fi


TMP_SNAKEFILE=$(python3 scripts/config_for_k8s.py -s Snakefile -c $CONFIG)


snakemake --kubernetes \
 --snakefile "$TMP_SNAKEFILE" \
 --default-remote-provider 'GS' \
 --default-remote-prefix "$BUCKET" \
 "$ARGS" && rm $TMP_SNAKEFILE || rm $TMP_SNAKEFILE

