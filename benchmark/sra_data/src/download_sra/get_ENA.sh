#!/bin/bash
set -euo pipefail

# usage: ./get_ENA.sh <ENA_report_file> <outdir>

ENA=$1
OUT=$(realpath $2)

if [ "$#" -eq 3 ]; then
    SSH="ssh $3"
else
		SSH=""
fi

mkdir -p ${OUT}

cut -f8 $ENA | awk -F';' -vOFS='\t' 'BEGIN {getline;}{print $1"\n"$2}' | parallel --jobs 1 $SSH bash ${PWD}/src/download_sra/get_one_ENA.sh {} $OUT

awk -F'\t' -vOFS='\t' 'BEGIN {getline} {split($7, a, ";"); split($8, b, ";");print a[1], b[1]
"\n" a[2], b[2]}' filereport_read_run_PRJNA606282_tsv.txt

python3 data/metadata/PRJNA606282/checksums.py $ENA > ${OUT}/checksums.md5


cd $OUT
md5sum -c checksums.md5

