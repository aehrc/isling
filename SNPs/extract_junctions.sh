#!/bin/bash


#inputs
allInt=$1 	#integration site txt file
hSing=$2	#human mapped bam file
vSing=$3	#virus mapped bam file

#intermediate files
hqueries=$4

#outputs
hSingJunctSam=$5	#human output sam
hSingJunctBam=$6	#human output bam
vSingJunctSam=$7	#virus output sam
vSingJunctBam=$8	#virus output bam

#integration site txt file format:
#field1: host chr
#field2: host int start
#field3: host int stop
#field4: virus ref
#field5: virus start
#field6: virus stop
#field19: read name

#get host coordinates and read names
cut -f1,2,3,4,5,6,19 $allInt | tail -n +2 > $hqueries

#if we had some integration sites
if [ -s $hqueries ]; then
	#get header from input bam
	samtools view -H $hSing > $hSingJunctSam
	samtools view -H $vSing > $vSingJunctSam
	#get reads from host 
	while IFS= read -r line; do
		hchr=`echo $line | awk '{print $1}'`
		hstart=`echo $line | awk '{print $2}'`
		hstop=`echo $line | awk '{print $3}'`
		vchr=`echo $line | awk '{print $4}'`
		vstart=`echo $line | awk '{print $5}'`
		vstop=`echo $line | awk '{print $6}'`
		name=`echo $line | awk '{print $7}'`
		samtools view $hSing chr${hchr}:${hstart}-${hstop} | grep $name >> $hSingJunctSam
		samtools view $vSing ${vchr}:${vstart}-${vstop} | grep $name >> $vSingJunctSam
	done < $hqueries
	samtools sort -o $hSingJunctBam $hSingJunctSam
	samtools sort -o $vSingJunctBam $vSingJunctSam
	samtools index $hSingJunctBam
	samtools index $vSingJunctBam
else
	touch $hSingJunctBam
	touch $hSingJunctSam
	touch $vSingJunctBam
	touch $vSingJunctSam
fi
