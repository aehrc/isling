#!/bin/bash

CORES=$1

mkdir -p data/reads/SRP023539

cd data/reads/SRP023539

ACCS="../../metadata/SRP023539/accs.txt"

cat $ACCS | parallel -j $CORES wget {}

bunzip2 *.bz2
