#!/bin/bash


## snakemake inputs
#input:
#	allInt = "../../out/{dset}/ints/{samp}.{host}.{virus}.integrations.txt",
#	vSing = "../../out/{dset}/virus_aligned_bam/{samp}.{virus}.bwa.bam"
#output:
#	vSingJunctSam = temp("../../out/{dset}/junct_aligns/{samp}.{virus}.bwa.junctReads.sam"),
#	vSingJunctBam = "../../out/{dset}/junct_aligns/{samp}.{virus}.bwa.junctReads.bam",
#	queries = temp("../../out/{dset}/ints/{samp}.{virus}.bwa.junctReads.txt")

##test files
#inputs
allInt=/scratch1/sco305/intvi_cmri/out/london_macaque_OTC/ints/7_S1.macaque.pAAV2-OTC.integrations.txt
vSing=/scratch1/sco305/intvi_cmri/out/london_macaque_OTC/virus_aligned_bam7_S1.pAAV2-OTC.bwa.mapped.bam
hSing=/scratch1/sco305/intvi_cmri/out/london_macaque_OTC/host_aligned_bam/7_S1.macaque.bwa.pAAV2-OTCMappedreads.bam

#intermediate
vqueries=7.vJunct.txt
hqueries=7.hJunct.txt

#outputs
vSingJunctSam=7.vMapped.junct.sam
vSingJunctBam=7.vMapped.junct.bam
hSingJunctSam=7.hMapped.junct.sam
hSingJunctBam=7.hMapped.junct.bam

#integration site txt file format:
#field1: host chr
#field2: host int start
#field3: host int stop
#field4: virus ref
#field5: virus start
#field6: virus stop
#field19: read name

#get host coordinates and read names
cut -f1,2,3,19 $allInt | tail -n +2 > $hqueries

#get virus coordinates and read names
cut -f4,5,6,19 $allInt | tail -n +2 > $vqueries


#if we had some integration sites
if [ -s $hqueries ]; then
	#get header from input bam
	samtools view -H $hSing > $hSingJunctSam
	#get reads from host 
	while IFS= read -r line; do
		chr=`echo $line | awk '{print $1}'`
		start=`echo $line | awk '{print $2}'`
		stop=`echo $line | awk '{print $3}'`
		name=`echo $line | awk '{print $4}'`
		samtools view $hSing chr${chr}:${start}-${stop} | grep $name >> $hSingJunctSam
	done < $hqueries
	samtools sort -o $hSingJunctBam $hSingJunctSam
	samtools index $hSingJunctBam
else
	touch $hSingJunctBam
	touch $hSingJunctSam
fi
