#!/bin/bash

#count mapped and unmapped reads using samtools


PROJ="`pwd`/../.."  #NOTE NEED TO RUN FROM scripts DIRECTORY IN THE CORRECT PROJECT!
DATA=${PROJ}/out
OUT=${PROJ}/out/summary/short_test/count_mapped.txt

printf "sample\tdataset\tmapped\tunmapped\n" > ${OUT}
files=($(find ${DATA}/*/*/*bam -maxdepth 0 -type f))
for file in "${files[@]}"; do
	DATASET=`dirname $file | xargs dirname | xargs basename`
	SAMPLE=`basename $file`
	MAPPED=`samtools view -c -F 4 -F 2048 $file` #not unmapped (-F 4), not supplementary (-F 2048)
	UNMAPPED=`samtools view -c -f 4 $file` #is unmapped (-f 4)
	printf "%s\t%s\t%s\t%s\n" $SAMPLE $DATASET $MAPPED $UNMAPPED >> ${OUT}
done



