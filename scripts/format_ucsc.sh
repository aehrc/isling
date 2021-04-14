#!/bin/bash

cd $1

DSET=$(pwd | xargs basename)
OUTFILE="${DSET}.post.bed"
#write header to each file containing name of sample and datatset
if [ -e "empty.bed" ]; then
	mv "empty.bed" "${OUTFILE}"
else
	FILES=*.bed
	for f in $FILES
	do
  		#add header to each set of data
  		echo "adding header to $f"
  		echo $(ls $f | perl -MEnv=DSET -ne '$_ =~ /\Q$DSET\E\.(.+)\.post\.bed/; print "track name=$1 description=\"integration sites from sample $1\"\n"') | cat - $f > "${f}.tmp"
  		rm $f
	done

	#combine samples for each dataset into one bed file
	echo "combining samples into ${OUTFILE}"
	
	cat *.bed.tmp > "${OUTFILE}"

	#clean up temp files
	rm *.tmp
	
fi

mv $OUTFILE ..
DIR=$(pwd | xargs basename)
cd .. && rmdir $DIR
