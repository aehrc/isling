#!/bin/bash

cd $1

#write header to each file containing name of sample and datatset
FILES=*.bed
for f in $FILES
do
  	#add header to each set of data
  	TEMP=$(echo ${f}.tmp)
  	echo `ls $f |  perl -ne '/(.+)(?=\.post\.bed)/; my @spl = split(/\./, $&); print "track name=@spl[1] description=\"integration sites from sample @spl[1], dataset @spl[0]\" \n" '` | cat - $f > $TEMP
  	rm $f
done

#combine samples for each dataset into one bed file
ls *.bed.tmp |  perl -ne '/(.+)(?=\..+\.post\.bed)/; print "$&\n"' | sort | uniq | parallel 'cat {}*.bed.tmp > {}.post.bed'

#clean up temp files
rm *.tmp
