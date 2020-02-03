#!/usr/bin/python

import sys, getopt
import os
import argparse
import subprocess
import re

def build_bwa_args(path, index, read1, read2, output, threshold, match, mismatch, h, threads, insert_open, insert_extend, delete_open, delete_extend):
	#build bwa args
	bwa_args = [path, "mem", "-t", str(threads), "-A", str(match), "-B", str(mismatch), "-O", "{},{}".format(delete_open, insert_open), "-E", "{},{}".format(delete_extend, insert_extend), "-L", "0,0", "-T", str(threshold), "-h", str(h), "-M", index, read1]

	#check if read1 and 2, or just read1
	if read2 is not None:
		bwa_args.append(read2)
	return bwa_args


def align_seqs(path, index, read1, read2, output, threshold, match, mismatch, h, threads, insert_open, insert_extend, delete_open, delete_extend):

	#if insert is a comma-seperated string of values, split into list
	if re.search(",", insert_open):
		insert_open = insert_open.split(",")
	
	#construct a list of bwa arguments, and a list of output filenames
	bwa_args_lst = []
	outputs = []

	#get part of output filename before ".sam"
	#first strip out .insert\d+ if it's in there
	output_name = re.search("(.+)\.insertOpen\d+(.+)", output)
	if output_name:
		#if there is a match, splice together two parts either side of insert
		output_name = output_name.group(1)+output_name.group(2)
	else:
		#otherwise just use output as provided to script
		output_name = output
	#remove .sam
	base_name = re.search("^(.+)*\.sam$", output_name)
	if base_name:
		output_base = base_name.group(1)
	else: 
		raise ValueError("Check output name contains .sam")

	#if multiple alingments to be done
	if len(insert) > 1:
		for val in insert:
			#add bwa args to list
			bwa_args_lst.append(build_bwa_args(path, index, read1, read2, output, threshold, match, mismatch, h, threads, val, insert_extend, delete_open, delete_extend))
			#add output to list
			outputs.append("{}.insert{}.sam".format(output_base, val))

	#if only one alignment
	else:
		bwa_args_lst.append(build_bwa_args(path, index, read1, read2, output, threshold, match, mismatch, h, threads, insert_open, insert_extend, delete_open, delete_extend))
		outputs.append(output)


	#do aligning 
	for i, args in enumerate(bwa_args_lst):
		p = subprocess.Popen(args, stdout=open(outputs[i], "w"))
		p.communicate()
		p.wait()

def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='Align single reads with BWA')
	parser.add_argument('--index', help='Index file for alignment',  required=True)
	parser.add_argument('--read1', help='Read 1 file', required=True)
	parser.add_argument('--read2', help='Read 2 file', required=False)
	parser.add_argument('--threshold', help='Alignment threshold', required=False, default='30')
	parser.add_argument('--match', help='Score for sequence match', required=False, default=1)
	parser.add_argument('--mismatch', help='Penalty for mismatch', required=False, default=2)
	parser.add_argument('--insert_open', help='Penalty for insertion open', required=False, default='6')
	parser.add_argument('--insert_extend', help='Penalty for insertion open', required=False, default='1')
	parser.add_argument('--delete_open', help='Pentalty for deletion open', required=False, default='6')
	parser.add_argument('--delete_extend', help='Penalty for insertion open', required=False, default='1')
	parser.add_argument('--bwa', help='Path to BWA', required=False, default='bwa')
	parser.add_argument('--output', help='Output', required=True)
	parser.add_argument('--hflag', help='-h value (INT): if there are <INT hits with score >80% of the max score, output all in XA', required=False, default="5")
	parser.add_argument('--threads', help='number of threads', required=False, default=1)
	args = parser.parse_args()

	#call to align_seqs
	align_seqs(args.bwa, args.index, args.read1, args.read2, args.output, args.threshold, args.match, args.mismatch, args.hflag, args.threads, args.insert_open, args.insert_extend, args.delete_open, args.delete_extend)

	

if __name__ == "__main__":
	main(sys.argv[1:])
