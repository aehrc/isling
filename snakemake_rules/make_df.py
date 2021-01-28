#### import modules ####
from glob import glob
from os import path, getcwd
import pandas as pd
import collections
import itertools
import pdb
import os
import subprocess
import scripts.pybam3 as pybam


#### defaults ####
bwa_mem_default = "-A 1 -B 2 -O 6,6 -E 1,1 -L 0,0 -T 10 -h 200"
merge_dist_default = 100
tol_default = 3
cutoff_default = 20
min_mapq_default = 0
merge_n_min_default = 1
split_default = 1
mean_frag_len_default = 'estimate'
align_cpus_default = 5

#### check file extensions ####

# supported file extensions for fastq files are .fastq and .fq, and compression with .gz and .bz2
# supported file extensions for bam/sam files are .bam, .sam

fastq_extensions = [".fq", ".fastq"]
bwa_prefix_extensions = [".amb", ".ann", ".bwt", ".pac", ".sa"]
sam_extensions = [".sam", ".bam"]
compression_extensions = [".gz", ".bz2"]

#### get information for wildcards ####

# get the name of each sample in each dataset, and save information about 
# how to process it in a dataframe

def make_df(config):

	# global options are specified with as a 'dataset' with the name 'global'
	# values in this dataset are applied to keys which are unset in the remaining datasets
	if 'global' in config:		
		# get default (global) options
		default = config.pop('global')
		for dataset in config:
			for key in default:
				if key not in config[dataset]:
					config[dataset][key] = default[key]


	rows = []

	for dataset in config:
	
		cat = check_read_suffixes(config, dataset)
		
		# get output directory
		outdir = get_value_or_default(config, dataset, 'out_dir', getcwd())
		outdir = path.normpath(outdir)
		config[dataset]['out_dir'] = outdir
		
		# get read directory
		readdir = check_required(config, dataset, 'read_folder')
		readdir = path.normpath(readdir)

		# Create temporary symlink directory
#		readdir = path.normpath(symlink_reads(readdir, os.path.join(outdir,"symreads"), config[dataset]['R1_suffix'], config[dataset]['R2_suffix']))
#		config[dataset]['read_folder'] = readdir
#		config[dataset]["R1_suffix"] = "_1.fq"
#		config[dataset]["R2_suffix"] = "_2.fq"

		# figure out if 'dedup', 'merge' and 'trim' are true or false for this dataset`
		dedup = check_bools(config, dataset, 'dedup')
		merge = check_bools(config, dataset, 'merge')
		if merge == 1:
			trim = 1
		else:
			trim = check_bools(config, dataset, 'trim')
		
		# get host and virus 
		# host and virus can either be spcified as 'host_name' (mandatory) and 
		# either 'host_prefix' or 'host_fasta' 
		# ('virus_name', 'virus_fasta', 'virus_prefix' in the case of viral reference), 
		
		# or in 'host'/'virus', where the key is the name of the host/virus and the 
		# value is the path the the fasta (for fasta input), or 'host_prefixes'/'virus_preifxes'
		# in the case of pre-indexed references
		
		get_refs(config, dataset, 'host')
		get_refs(config, dataset, 'virus')		
		
		# only need adapters if we're trimming or merging
		if merge or trim:
			adapter_1 =check_required(config, dataset, 'read1-adapt')
			adapter_2 =check_required(config, dataset, 'read2-adapt')
		else:
			adapter_1 = ""
			adapter_2 = ""
		
		# check for other optional parameters
		mean_frag_len = get_value_or_default(config, dataset, 'mean-frag-len', mean_frag_len_default)
		merge_dist = get_value_or_default(config, dataset, 'merge-dist', merge_dist_default)
		clip_cutoff = get_value_or_default(config, dataset, 'clip-cutoff', cutoff_default)
		cigar_tol = get_value_or_default(config, dataset, 'cigar-tol', tol_default)
		min_mapq = get_value_or_default(config, dataset, 'min-mapq', min_mapq_default)
		bwa_mem_params = get_value_or_default(config, dataset, 'bwa-mem', bwa_mem_default)
		merge_n_min = get_value_or_default(config, dataset, 'merge-n-min', merge_n_min_default)
		split = get_value_or_default(config, dataset, 'split', split_default)
		align_cpus = get_value_or_default(config, dataset, 'align-cpus', align_cpus_default)
		
		# check values are integers greater than a minimum value
		check_int_gt(merge_dist, -1, 'merge-dist', dataset)
		check_int_gt(clip_cutoff, 1, 'clip-cutoff', dataset)
		check_int_gt(cigar_tol, 0, 'cigar-tol', dataset)
		check_int_gt(min_mapq, -1, 'min-mapq', dataset)
		check_int_gt(merge_n_min, -1, 'merge-n-min', dataset)
		check_int_gt(split, 0, 'split', dataset)

		check_int_gt(min_mapq, -1, 'merge-n-min', dataset)
		if mean_frag_len != 'estimate':
			if mean_frag_len < 0:
				raise ValueError('key "mean-frag-len" in dataset {dataset} should be either the string "estimate", or a number greater than 0')
		elif isinstance(mean_frag_len, str):
			if mean_frag_len != 'estimate':
				raise ValueError('key "mean-frag-len" in dataset {dataset} should be either the string "estimate", or a number greater than 0')
		
		# get arguments for running postprocessing scripts
		postargs = make_post_args({dataset : config[dataset]})[0][dataset]
		
		# get samples for this dataset
		samples, is_bam = get_samples(config, dataset)
		
		# make one row for each sample
		for sample, host, virus in itertools.product(samples,
																									config[dataset]['host'].keys(), 
																									config[dataset]['virus'].keys()):
			print(is_bam)
			if is_bam:
				bam_file = f"{readdir}/{sample}{config[dataset]['bam_suffix']}"
				R1_file = f"{outdir}/{dataset}/reads/{sample}{config[dataset]['R1_suffix']}"
				R2_file = f"{outdir}/{dataset}/reads/{sample}{config[dataset]['R2_suffix']}"
				print(R2_file)
				split_lines = split_lines_bam(bam_file, split, cat)
			else:
				bam_file = ""	
				R1_file = f"{readdir}/{sample}{config[dataset]['R1_suffix']}"
				R2_file = f"{readdir}/{sample}{config[dataset]['R2_suffix']}"
				split_lines = split_lines_fastq(R1_file, split, cat)
			
			# if there is more than one virus or host, append these to the dataset name
			if len(config[dataset]['host'].keys()) > 1 or len(config[dataset]['virus'].keys()) > 1:
				dataset_name = f"{dataset}_{host}_{virus}"
			else:
				dataset_name = dataset
				
			# make sample-specific information
			unique = f"{dataset_name}+++{sample}"
			
			host_fasta = config[dataset]["host"][host]
			host_prefix = config[dataset]["host_prefixes"][host]
			virus_fasta = config[dataset]["virus"][virus]
			virus_prefix = config[dataset]["virus_prefixes"][virus]	
			
			# append combinations of each sample, host and virus		
			rows.append((dataset_name, dataset, sample, host, host_fasta, host_prefix, virus, virus_fasta, virus_prefix, merge, trim, dedup, unique,  outdir, bwa_mem_params, R1_file, R2_file, bam_file, adapter_1, adapter_2, postargs, merge_dist, merge_n_min, clip_cutoff, cigar_tol, min_mapq, split, mean_frag_len, align_cpus, cat, split_lines))

			
	# check there aren't any duplicate rows
	if len(set(rows)) != len(rows):
		raise ValueError("Error - configfile results in duplicate analyses, check samples and dataset names are unique")
	
	# make dataframe
	toDo = pd.DataFrame(rows, columns=['dataset', 'config_dataset', 'sample', 'host', 'host_fasta', 'host_prefix', 'virus', 'virus_fasta', 'virus_prefix', 'merge', 'trim', 'dedup', 'unique', 'outdir', 'bwa_mem_params', 'R1_file', 'R2_file', 'bam_file', 'adapter_1', 'adapter_2', 'postargs', 'merge_dist', 'merge_n_min', 'clip_cutoff', 'cigar_tol', 'min_mapq', 'split', 'mean_frag_len', 'align_cpus', 'cat', 'split_lines'])
	
	# do checks on dataframe
	check_dataset_sample_unique(toDo)
	
	ref_names = make_reference_dict(toDo)
	check_fastas_unique(toDo, ref_names)

	print(toDo)
	return toDo


def check_dataset_sample_unique(toDo):
	# check that every combination of 'dataset' and 'sample' is unique
	if len(set(toDo.loc[:,'unique'])) != len(toDo.loc[:,'unique']):
		raise ValueError("Every combination of 'dataset' and 'sample' must be unique! Check inputs and try again")

def check_fastas_unique(toDo, ref_names):
	# check that each host/virus name refers to only one fasta
	for i, virus_name in enumerate(toDo.loc[:,'virus']):
		if toDo.loc[i,'virus_fasta'] != "":
			if toDo.loc[i,'virus_fasta'] != ref_names[virus_name]:
				raise ValueError(f"Virus {virus_name} is used as a name in multiple datasets for different fasta files.  Please modify your configfile so that each unique virus/host name only refers to one fasta file")
	for i, host_name in enumerate(toDo.loc[:,'host']):
		if toDo.loc[i,'host_fasta'] != "":
			if toDo.loc[i,'host_fasta'] != ref_names[host_name]:
				raise ValueError(f"Host {host_name} is used as a name in multiple datasets for different fasta files.  Please modify your configfile so that each unique virus/host name only refers to one fasta file")
	
	
def make_reference_dict(toDo):
	# construct dictionary with reference names as keys and reference fastas as values
	
	ref_names = {}
	for index, row in toDo.iterrows():
		if row['host_fasta'] != "":
			ref_names[row['host']] = row['host_fasta']
		if row['virus_fasta'] != "":
			ref_names[row['virus']] = row['virus_fasta']
	
	return ref_names
	
# construct arguments for postprocess.R script for each dataset
def make_post_args(config):

	POSTARGS = {}
	TOSORT = []
	SORTED = []
	for dataset in config:
		POSTARGS[dataset] = []
		if "post" in config[dataset]:
			for element in config[dataset]["post"]:
				# need to check if this element is a string or a dict
				if isinstance(element, str):
					# look for keys to be 'filter' or 'dedup'
					if element == "filter":
						POSTARGS[dataset].append("filter")
					elif element == "dedup":
						POSTARGS[dataset].append("dedup")
				# postprocessing types with files specified will be in ordered dictionaries
				elif isinstance(element, collections.OrderedDict):
					if "mask-exclude" in element.keys():
						for bed in element["mask-exclude"]:
							POSTARGS[dataset].append("mask-exclude")
							POSTARGS[dataset].append(bed)	
					elif "mask-include" in element.keys():
						for bed in element["mask-include"]:
							POSTARGS[dataset].append("mask-include")
							POSTARGS[dataset].append(bed)
					elif "nearest-gtf" in element.keys():
						for gtf in element["nearest-gtf"]:
							sortedgtf = path.splitext(gtf)[0] + ".sorted.gtf"
							POSTARGS[dataset].append("nearest-gtf")
							POSTARGS[dataset].append(sortedgtf)
							if gtf not in TOSORT:
								TOSORT.append(gtf)
								SORTED.append(sortedgtf)
					elif "nearest-bed" in element.keys():
						for bed in element["nearest-bed"]:
							sortedbed = path.splitext(bed)[0] + ".sorted.bed"
							POSTARGS[dataset].append("nearest-bed")
							POSTARGS[dataset].append(sortedbed)
							if bed not in TOSORT:
								TOSORT.append(bed)
								SORTED.append(sortedbed)
					elif "RNA-seq" in element.keys():
						ref = element["genes"]
						sortedref = path.splitext(ref)[0] + ".sorted" + path.splitext(ref)[1]
						if ref not in TOSORT:
							TOSORT.append(ref)
							SORTED.append(sortedref)
						for tsv in element["counts"]:
							POSTARGS[dataset].append("RNA-seq-gtf")
							POSTARGS[dataset].append(sortedref)
							POSTARGS[dataset].append(element["col"])
							POSTARGS[dataset].append(tsv)
			POSTARGS[dataset] = " ".join(POSTARGS[dataset])
		
	return POSTARGS, TOSORT, SORTED
	

def get_refs(config, dataset, ref_type):
	"""
	host and virus can either be spcified as 'host_name' (mandatory) and 
	either 'host_prefix' or 'host_fasta' 
	('virus_name', 'virus_fasta', 'virus_prefix' in the case of viral reference), 
		
	or in 'host'/'virus', where the key is the name of the host/virus and the 
	value is the path the the fasta (for fasta input), or 'host_prefixes'/'virus_prefixes'
	in the case of pre-indexed references
	
	if both prefixes and fastas are specified, use only prefixes and ignore fastas
	
	make sure the necessary keys are present in config for this dataset
	"""
	assert ref_type in {'host', 'virus'}
	# if one prefix specified
	if f"{ref_type}_prefix" in config[dataset]:
		prefix = config[dataset][f"{ref_type}_prefix"]
		# host name is mandatory
		name = check_required(config, dataset, f"{ref_type}_name")
		
		# add dict to config
		config[dataset][f"{ref_type}_prefixes"] = {name:prefix}
		
		# add fasta info
		if f"{ref_type}_fasta" in config[dataset]:
			print(f"dataset {dataset} has both fasta and bwa prefix specified: ignoring fasta")
			
		config[dataset][ref_type] = {name:""}
		
		return

	# if prefixes specified in a dict
	elif f"{ref_type}_prefixes" in config[dataset]:
		# make sure dict contains at least one entry
		if len(config[dataset][f"{ref_type}_prefixes"]) == 0:
			raise ValueError(f"At least one {ref_type} prefix must be specified for dataset {dataset}")
		
		# check for fasta dict, print warning if exists, and create fasta dict with empty strings
		if ref_type in config[dataset]:
			print(f"dataset {dataset} has both fastas and bwa prefixes specified: ignoring fastas")			
		
		config[dataset][ref_type] = {name:"" for name in config[dataset][f"{ref_type}_prefixes"]}
		
	# if one fasta specified	
	elif f"{ref_type}_fasta" in config[dataset]:	
		fasta = config[dataset][f"{ref_type}_fasta"]
		# host name is mandatory
		name = check_required(config, dataset, f"{ref_type}_name")
		
		# add dict to config
		config[dataset][ref_type] = {name:fasta}	
		
		# add prefix info
		outdir = config[dataset]['out_dir']
		prefix = f"{outdir}/references/{name}/{name}"
		
		# add dict to config
		config[dataset][f"{ref_type}_prefixes"] = {name:prefix}	
		return	
		
	# fasta refs specified in a dict
	elif f"{ref_type}" in config[dataset]:
		# make sure dict contains at least one entry
		if len(config[dataset][f"{ref_type}"]) == 0:
			raise ValueError(f"At least one {ref_type} fasta must be specified for dataset {dataset}")
		
		outdir = config[dataset]['out_dir']
		config[dataset][f"{ref_type}_prefixes"] = {name:f"{outdir}/references/{name}/{name}" for name in config[dataset][f"{ref_type}"]}
		


def check_bools(config, dataset, key):
	"""
	Check that the key in the speicfied dataset is defined, and return 1 if true, 0 if false
	"""
	# if not specified
	if key not in config[dataset]:
		raise ValueError(f"Please specify True or False for {key} in dataset {dataset}")
	
	# try to figure out if the user wanted true or false
	if isinstance(config[dataset][key], bool):
		return int(config[dataset][key])
	elif isinstance(config[dataset][key], str):
		if config[dataset][key].lower() == "true":
			return 1
		else:
			return 0
	else:
		raise ValueError(f"Please specify True or False for '{key}' in dataset {dataset}")
	
def check_required(config, dataset, key):
	"""
	Check that a required key has been specified, and raise an error if it hasn't
	"""
	if key not in config[dataset]:
		raise ValueError(f"Please specify required parameter '{key}' in dataset {dataset}")
		
	return config[dataset][key]
		
def get_value_or_default(config, dataset, key, default):
	"""
	Check that a optional paramter has been specified, and use the default if it hasn't
	"""
	
	if key not in config[dataset]:
		print(f"Paramter '{key}' not specified for dataset {dataset}: using default {default}")
		return default
	
	return config[dataset][key]	
	
def check_int_gt(var, minimum, name, dataset):
	"""
	check if a variable is an integer greater than 'minimum'
	"""	
	if not isinstance(var, int):
		raise ValueError(f"Parameter {name} in dataset {dataset} must be an integer")
	
	if var < minimum:
		raise ValueError(f"Parameter {name} in dataset {dataset} must be greater than {minimum}")

def check_read_suffixes(config, dataset):

	cat = 'cat'

	# if input is fastq
	if "R1_suffix" in config[dataset]:
	
		# check that R2_suffix is also specified
		if "R2_suffix" not in config[dataset]:
			raise ValueError("If R1_suffix is specified (for fastq input), R2_suffix must also be specified")
				
		if config[dataset]['R1_suffix'] == config[dataset]['R2_suffix']:
			raise ValueError(f"R1_suffix must be different to R2_suffix for dataset {dataset}")
	
		# check R1_suffix for compression and fastq file extension
		R1_first_extension = path.splitext(config[dataset]["R1_suffix"])[1]
		R2_first_extension = path.splitext(config[dataset]["R2_suffix"])[1]
		
		#bz2 or gz compression
		if R1_first_extension in compression_extensions:
		
			# check for which cat we should use
			cat =  'zcat' if R1_first_extension == '.gz' else 'bzcat'
		
			# get second extensions
			R1_second_extension = path.splitext(path.splitext(config[dataset]["R1_suffix"])[0])[1]
			R2_second_extension = path.splitext(path.splitext(config[dataset]["R2_suffix"])[0])[1]
		
			# check that input looks like a fastq files
			if R1_second_extension not in fastq_extensions or R2_second_extension not in fastq_extensions:
				raise InputError("input files do not look like fastq files (extension is not '.fq' or '.fastq'")

			
		# if uncompressed, check file looks like fastq file
		elif R1_first_extension not in fastq_extensions or R2_first_extension not in fastq_extensions:
			raise InputError("for fastq files, input file extensions must be '.fastq' or '.fq'")

	# if input is bam
	elif "bam_suffix" in config[dataset]:
		extension = config[dataset]["bam_suffix"]
		if extension != ".bam" and extension != ".sam":
			extension = path.splitext(config[dataset]["bam_suffix"])[1]
			if extension != ".bam" and extension != ".sam":
				raise InputError("For aligned input, only '.bam' and '.sam' files are currently supported")
	
	# if nether R1/R2 suffix or bam suffix specified
	else:
		raise InputError("Please specify either 'R1_suffix' and 'R2_suffix' or 'bam_suffix' in the config file")


	return cat

def get_samples(config, dataset):
		
		# if samples are specified
		if 'samples' in config[dataset]:
			if not hasattr(config[dataset]['samples'], '__iter__'):
				raise ValueError(f"the value for key 'samples' in dataset {dataset} should be a list")

			if 'R1_suffix' in config[dataset]:
				is_bam = False
			elif 'bam_suffix' in config[dataset]:
				is_bam = True
				config[dataset]["R1_suffix"] = "_1.fq.gz"
				config[dataset]["R2_suffix"] = "_2.fq.gz"
			else:
				raise ValueError(f"please specfiy either 'R1_suffix' and 'R2_suffix' or 'bam_suffix' for dataset {dataset}")
				
			return config[dataset]['samples'], is_bam
			
		
		# get fastq files for input
		if "R1_suffix" in config[dataset]:
			
			is_bam = False
			
			# get samples with names ending in R1_suffix
			suffix = config[dataset]["R1_suffix"]
			config[dataset]['read_folder'] = path.normpath(config[dataset]['read_folder'])
			folder = config[dataset]['read_folder']	
			samples = [path.basename(f)[:-len(suffix)] for f in glob(f"{folder}/*{suffix}")]
			
			# check for corresponding R2 files
			dropped_samples = []
			R2_suffix = config[dataset]["R2_suffix"]
			for R2_file in [f"{folder}/{sample}{R2_suffix}" for sample in samples]:
				if not path.exists(R2_file):
					print(f"Found R1 file for sample {sample} in dataset {dataset}, bud didn't find matching R2 file.  Ignoring R1 file")
					dropped_samples.append(sample)
			# drop samples that we couldn't find matching R1 and R2 files for
			samples = [samp for samp in samples if samp not in dropped_samples]
			
			if len(samples) == 0:
				print(f"warning: no files found for dataset {dataset}")
		
		# get bam/sam files for input
		elif "bam_suffix" in config[dataset]:
			is_bam = True
			suffix = config[dataset]["bam_suffix"]
			folder = config[dataset]['read_folder']	
			config[dataset]["R1_suffix"] = "_1.fq.gz"
			config[dataset]["R2_suffix"] = "_2.fq.gz"
			samples = [path.basename(f)[:-len(suffix)] for f in glob(f"{folder}/*{suffix}")]
			if len(samples) == 0:
				print(f"warning: no files found for dataset {dataset}")
		else:
			raise InputError(f"please specify either 'bam_suffix' or 'R1_suffix' and 'R2_suffix' for dataset {dataset}")
			
		return samples, is_bam


def split_lines_fastq(read_file_path, split_n, cat):
	# in the trivial case, only one part
	if split_n == 1:
		return ''

	# get number of lines in file
	if cat == 'bzcat':
		count = int(subprocess.check_output(['bzgrep', '-Ec', '$', read_file_path]).split()[0])
	elif cat == 'zcat':
		count = int(subprocess.check_output('zgrep', '-Ec', '$', read_file_path]).split()[0])
	else:
		count = int(subprocess.check_output(['wc', '-l', read_file_path]).split()[0])	
	print(f"Counted {count} lines in {read_file_path}")

	# we want to split the line into split_n chunks
	# where each chunk starts on a line (1+4i)
	# and each chunck ends on a line (4j)
	# where i, j are integers
	chunk_lines = count / split_n
	chunks = []
	curr_start = 1
	while curr_start < count:
		end = 4 * round((curr_start + chunk_lines) / 4)
		chunks.append(f"{curr_start},{end}")
		curr_start = end + 1
		
	return "xxx".join(chunks)


def split_lines_bam(bam_file_path, split_n, cat):
	# in the trivial case, only one part
	if split_n == 1:
		return ''

	# get number of lines in file
	alignments = pybam.read(bam_file_path)
	count = 0
	for alignment in alignments:
		count += 1

	header_count = alignments.file_header.split('\n')
	header_count = list(filter(None, header_count))

	# we want to split the line into split_n chunks
	# where each chunk starts on a line (1+4i)
	# and each chunck ends on a line (4j)
	# where i, j are integers
	chunk_lines = count / split_n
	chunks = []
	curr_start = 1
	while curr_start < count:
		end = 4 * round((curr_start + chunk_lines) / 4)
		chunks.append(f"{curr_start},{end}")
		curr_start = end + 1
	
	return "xxx".join(chunks)

