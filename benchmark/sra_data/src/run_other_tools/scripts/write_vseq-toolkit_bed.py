#!/usr/bin/env python3

import argparse
import sys
import csv
import pdb

def main(args):

	parser = argparse.ArgumentParser(description = "convert seeksv output to bed file")
	parser.add_argument('--vseq-toolkit-output', help = 'output from VSeq-Toolkit', required=True)
	parser.add_argument('--output', help = 'name of output bed file', default='vseq_toolkit.bed')
	args = parser.parse_args()

	# read lines of output and write bed file
	with open(args.vseq_toolkit_output, 'r', newline='') as vseq, open(args.output, 'w', newline='') as bed:
		reader = csv.DictReader(vseq, delimiter=',')
		writer = csv.writer(bed, delimiter='\t')
		
		for row in reader:
	
			if row['GenomicPosition'] == '':
				continue
			# https://github.com/CompMeth/VSeq-Toolkit/issues/4
			# OverlapFusion - the analysis strategy depends on finding a fusion between vector and genome mapped regions of the read. This column indicates how many bases at the fusion site were overlapped between genomic and vector mapped regions.
			#DistanceFusion - similarly it indicates if there were any bases at the fusion site between vector and genome that were not mapped.
			start = int(row['GenomicPosition'])
			stop = int(row['GenomicPosition']) + int(row['OverlapFusion']) + int(row['DistanceFusion'])
			
			#write_row = (row['Chr'], start, stop, row['StrandGenomic'])
			write_row = (row['Chr'], start, stop)
			writer.writerow(write_row)
	
	
if __name__ == "__main__":
	main(sys.argv)
