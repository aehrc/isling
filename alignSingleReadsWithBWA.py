#!/usr/bin/python

import sys, getopt
import os
import argparse
import subprocess

def align_seqs(path, index, reads, output, threshold, match, mismatch, aflag, h):	
	if (aflag =="True"):
		bwa_args = [path, "mem", "-A", str(match), "-B", str(mismatch), "-O", "6,6", "-L", "0,0", "-T", str(threshold), "-h", str(h), "-a", index, reads]
	else:
		bwa_args = [path, "mem", "-A", str(match), "-B", str(mismatch), "-O", "6,6", "-L", "0,0", "-T", str(threshold), "-h", str(h), index, reads]
	
	p = subprocess.Popen(bwa_args, stdout=open(output, 'w'))
	(output, err) = p.communicate()

def main(argv):
	parser = argparse.ArgumentParser(description='Align single reads with BWA')
	parser.add_argument('--index', help='Index file for alignment',  required=True)
	parser.add_argument('--reads', help='Read file', required=True)
	parser.add_argument('--threshold', help='Alignment threshold', required=False, default='30')
	parser.add_argument('--match', help='Score for sequence match', required=False, default=1)
	parser.add_argument('--mismatch', help='Penalty for mismatch', required=False, default=2)
	parser.add_argument('--bwa', help='Path to BWA', required=False, default='bwa')
	parser.add_argument('--output', help='Output', required=True)
	parser.add_argument('--aflag', help='-a flag set?', required=False, default="False", choices=["False", "True"])
	parser.add_argument('--hflag', help='-h value (INT): if there are <INT hits with score >80% of the max score, output all in XA', required=False, default="5")
	args = parser.parse_args()

	align_seqs(args.bwa, args.index, args.reads, args.output, args.threshold, args.match, args.mismatch, args.aflag, args.hflag)

if __name__ == "__main__":
	main(sys.argv[1:])
