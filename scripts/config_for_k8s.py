#!/usr/bin/env python3

## this script reads a yaml file, and adds it to the top of the isling snakefile

from sys import argv, exit
import os
from tempfile import NamedTemporaryFile
import argparse
import yaml
import pdb

def main(args):

	# parse arguments
	parser = argparse.ArgumentParser(description = "This script reads a config"
		" file and adds it to the top of a snakefile")
	parser.add_argument('--config', '-c', help = "Config file to read", 
							required=True)
	parser.add_argument('--snakefile', '-s', help = "Snakefile to add config to",
							required=True)
	args = parser.parse_args(args[1:])
	
	# check that input files exist
	check_file_exists(args.config)
	check_file_exists(args.snakefile)
	
	# read yaml
	with open(args.config, 'r') as stream:
		try:
			config = yaml.safe_load(stream)
		except yaml.YAMLError as exc:
			print("invalid config")
			print(exc)
			exit(1)
			
	if 'snakedir' in config:
		del config['snakedir']

	# get path of Snakefile for output of temporary file
	out_dir = os.path.dirname(os.path.realpath(args.snakefile))
	
	# make temporary file for writing snakefile
	with NamedTemporaryFile(mode = 'w', delete=False, dir = out_dir) as f:	
		
		# write config
		f.write(f"config = {str(config)}\n")
		
		# write rest of snakefile
		with open(args.snakefile, 'r') as snakefile:
			for line in snakefile:
				f.write(line)
		
		# print name of temp file with config
		print(os.path.basename(f.name))
	
			
def check_file_exists(file):
	if not os.path.exists(file):
		print(f"{file} does not exist")
		exit(1)	

if __name__ == "__main__":
	main(argv)
