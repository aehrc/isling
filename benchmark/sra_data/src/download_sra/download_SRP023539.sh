#!/bin/bash

CORES=$1

mkdir -p data/reads/SRP023539

cd data/reads/SRP023539

ACCS="../../data/metadata/SRP023539/accs.txt"

cat $ACCS | parallel -j $CORES wget {}

if [ ! -e SRR873836_2.fastq ] ; then
	bunzip2 *.bz2
fi
