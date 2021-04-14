#!/bin/bash
set -euo pipefail

# $1 - config file
# $2 - name of bucket
# other args - passed to snakemake (no need for --kubernetes or bucket details)

# script to run isling using kubernetes.  Pass through any other arguments passed the this script to snakemake

if [ $# -lt 2 ]; then
	echo "usage: ./run_kubernetes <config> <gs_bucket_name>"
	exit 1
fi

CONFIG=$1
BUCKET=$2

TMP_SNAKEFILE=$(python3 scripts/config_for_k8s.py -s Snakefile -c $CONFIG)
SMK="snakemake --kubernetes --snakefile $TMP_SNAKEFILE --default-remote-provider GS --default-remote-prefix $BUCKET"

if [ $# -gt 2 ]; then
	SMK="${SMK} ${@:3}"
fi

$SMK && rm $TMP_SNAKEFILE || rm $TMP_SNAKEFILE

