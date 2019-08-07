#!/bin/bash

#count number of potential short integrations
#pass as input path to output file
OUT=$1

echo "" > ../../out/summary/count_short.txt

for f in ../../out/*/ints/*short.txt; do
	echo "${f}	`tail -n+2 $f | wc -l`" >> ${OUT}
done
