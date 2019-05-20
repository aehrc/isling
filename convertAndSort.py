#!/usr/bin/python

import sys, getopt
import os
import argparse
import pysam

def main(argv):
	parser = argparse.ArgumentParser(description='convert and sort alignment files')
	parser.add_argument('--sam', help='sam file', required=True)
	parser.add_argument('--bam', help='bam file', required=True)
	parser.add_argument('--sort', help='sorted file', required=True)
	args = parser.parse_args()

	pysam.view("-bSo", args.bam, args.sam, catch_stdout=False)
	pysam.sort("-o", args.sort, args.bam)
	pysam.index(args.sort)

if __name__ == "__main__":
	main(sys.argv[1:])
