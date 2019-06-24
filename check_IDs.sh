#!/bin/bash

#check that all reads present in previous pipeline results are also present in new results
#pass in two directories with files to check as arguments

NEW=$1
OLD=$2

if [ ! $# == 2 ]; then
  echo "Usage: $0 <new_directory> <old_directory>"
  exit
fi

if [ ! -d ${NEW} ]; then
	echo "directory $NEW does not exist"
elif [ ! -d ${OLD} ]; then
	echo "directory $OLD does not exist"
fi
	

for filename in ${OLD}/*.txt; do
	basefile=$(basename ${filename})
	if [ ! -f ${OLD}/${basefile} ]; then
		echo "file ${filename} does not exist in new directory"
	else 
		echo "READ IDs differing for ${basefile}:"
		diff <(awk -F'\t' '{if (NR!=1) print $(NF-1)}' ${OLD}/${basefile} | sort ) <(awk -F'\t' '{if (NR!=1) print $(NF-1)}' ${NEW}/${basefile} | sort)
	fi
done
		

