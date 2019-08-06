#!/usr/bin/python

import sys, getopt
import os
import argparse
import subprocess

def align_seqs(path, index, read1, read2, output, threshold, match, mismatch, h, threads, insert, delete):

	#build bwa args
	bwa_args = [path, "mem", "-t", str(threads), "-A", str(match), "-B", str(mismatch), "-O", "{delete},{insert}", "-L", "0,0", "-T", str(threshold), "-h", str(h), index, read1]

	#check if read1 and 2, or just read1
	if read2 is not None:
		bwa_args.append(read2)


	#do aligning 
	p = subprocess.Popen(bwa_args, stdout=open(output, "w"))
	p.communicate()


def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='Align single reads with BWA')
	parser.add_argument('--index', help='Index file for alignment',  required=True)
	parser.add_argument('--read1', help='Read 1 file', required=True)
	parser.add_argument('--read2', help='Read 2 file', required=False)
	parser.add_argument('--threshold', help='Alignment threshold', required=False, default='30')
	parser.add_argument('--match', help='Score for sequence match', required=False, default=1)
	parser.add_argument('--mismatch', help='Penalty for mismatch', required=False, default=2)
	parser.add_argument('--insert', help='Penalty for insertion open', required=False, default=6)
	parser.add_argument('--delete', help='Pentalty for deletion open', required=False, default=6)
	parser.add_argument('--bwa', help='Path to BWA', required=False, default='bwa')
	parser.add_argument('--output', help='Output', required=True)
	parser.add_argument('--hflag', help='-h value (INT): if there are <INT hits with score >80% of the max score, output all in XA', required=False, default="5")
	parser.add_argument('--threads', help='number of threads', required=False, default=1)
	args = parser.parse_args()

	#call to align_seqs
	align_seqs(args.bwa, args.index, args.read1, args.read2, args.output, args.threshold, args.match, args.mismatch, args.hflag, args.threads, args.insert, args.delete)

	

if __name__ == "__main__":
	main(sys.argv[1:])
