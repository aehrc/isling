#### import modules ####
from glob import glob
from os import path, getcwd
import pandas as pd
import collections
import itertools
import pdb

#### custom errors ####

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message


#### defaults ####
bwa_mem_params="-A 1 -B 2 -O 6,6 -E 1,1 -L 0,0 -T 10 -h 200"


#### check file extensions ####

# supported file extensions for fastq files are .fastq and .fq, and compression with .gz and .bz2
# supported file extensions for bam/sam files are .bam, .sam

fastq_extensions = [".fq", ".fastq"]
sam_extensions = [".sam", ".bam"]
compression_extensions = [".gz", ".bz2"]

def check_input_files(config):

	for dataset in config:

		# if input is fastq
		if "R1_suffix" in config[dataset]:
	
			# check that R2_suffix is also specified
			if "R2_suffix" not in config[dataset]:
				raise InputError("If R1_suffix is specified (for fastq input), R2_suffix must also be specified")
	
			# check R1_suffix for compression and fastq file extension
			R1_first_extension = path.splitext(config[dataset]["R1_suffix"])[1]
			R2_first_extension = path.splitext(config[dataset]["R2_suffix"])[1]
		
			#bz2 or gz compression
			if R1_first_extension in compression_extensions:
		
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
		

#### get information for wildcards ####

# get the name of each sample in each dataset, and save information about 
# how to process it in a dataframe

def make_df(config):

	check_input_files(config)

	rows = []

	for dataset in config:
		# get output directory
		if "out_dir" not in config[dataset]:
			print(f"output directory not specified for dataset {dataset}, using current directory")
			config[dataset]["out_dir"] = getcwd()
		else:
			config[dataset]["out_dir"] = path.normpath(config[dataset]["out_dir"])
		
		# get fastq files for input
		if "R1_suffix" in config[dataset]:
			is_bam = False
			suffix = config[dataset]["R1_suffix"]
			config[dataset]['read_folder'] = path.normpath(config[dataset]['read_folder'])
			folder = config[dataset]['read_folder']	
			samples = [path.basename(f)[:-len(suffix)] for f in glob(f"{folder}/*{suffix}")]
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
			raise InputError(f"please specify either bam_suffix or R1_suffix and R2_suffix for dataset {dataset}")
			
		# figure out if 'dedup' and 'merge' are true or false for this dataset
		if isinstance(config[dataset]["dedup"], bool):
			dedup = int(config[dataset]["dedup"])
		elif isinstance(config[dataset]["dedup"], str):
			if config[dataset]["dedup"].lower() == "true":
				dedup = 1
			else:
				dedup = 0
		else:
			raise InputError(f"Please specify True or False for 'dedup' in dataset {dataset}")
		if isinstance(config[dataset]["merge"], bool):
			merge = int(config[dataset]["merge"])
		elif isinstance(config[dataset]["merge"], str):
			if config[dataset]["merge"].lower() == "true":
				merge = 1
			else:
				merge = 0
		else:
			raise InputError(f"Please specify True or False for 'merge' in dataset {dataset}")
		
		# get host and virus 
		# host and virus can either be spcified as 'host_name' and 'host_fasta' ('virus_name' and 'virus_fasta'), 
		# or in 'host'/'virus', where the key is the name of the host/virus and the value is the path the the fasta
		if 'host' in config[dataset]:
			
			# check at least one host was specified
			if len(config[dataset]['host']) == 0:
				raise ValueError(f"At least one host must be specified for dataset {dataset}")
			
		else:
			config[dataset]['host'] = {config[dataset]["host_name"] : config[dataset]["host_fasta"]}
			
		if 'virus' in config[dataset]:
			
			# check at least one virus was specified
			if len(config[dataset]['virus']) == 0:
				raise ValueError(f"At least one virus must be specified for dataset {dataset}")
			
		else:
			config[dataset]['virus'] = {config[dataset]["virus_name"] : config[dataset]["virus_fasta"]}
		
		adapter_1 = config[dataset]["read1-adapt"]
		adapter_2 = config[dataset]["read2-adapt"]
		
		# get arguments for running postprocessing scripts
		postargs = make_post_args({dataset : config[dataset]})[0][dataset]
		
		# make one row for each sample
		for sample, host, virus in itertools.product(samples,config[dataset]['host'].keys(), config[dataset]['virus']):
			
			if is_bam:
				bam_file = f"{path.normpath(config[dataset]['read_folder'])}/{sample}{config[dataset]['bam_suffix']}"
				R1_file = f"{path.normpath(config[dataset]['out_dir'])}/{dataset}/reads/{sample}{config[dataset]['R1_suffix']}"
				R2_file = f"{path.normpath(config[dataset]['out_dir'])}/{dataset}/reads/{sample}{config[dataset]['R2_suffix']}"
			else:
				bam_file = ""	
				R1_file = f"{path.normpath(config[dataset]['read_folder'])}/{sample}{config[dataset]['R1_suffix']}"
				R2_file = f"{path.normpath(config[dataset]['read_folder'])}/{sample}{config[dataset]['R2_suffix']}"
			
			# if there is more than one virus or host, append these to the dataset name
			if len(config[dataset]['host'].keys()) > 1 or len(config[dataset]['virus'].keys()) > 1:
				dataset_name = f"{dataset}_{host}_{virus}"
			else:
				dataset_name = dataset
				
			# make sample-specific information
			unique = f"{dataset_name}+++{sample}"
			
			# append combinations of each sample, host and virus		
			rows.append((dataset_name, dataset, sample, host, config[dataset]["host"][host], virus, config[dataset]["virus"][virus], merge, dedup, unique,  config[dataset]["out_dir"], bwa_mem_params, R1_file, R2_file, bam_file, adapter_1, adapter_2, postargs))

			
	# check there aren't any duplicate rows
	if len(set(rows)) != len(rows):
		raise ValueError("Error - configfile results in duplicate analyses, check samples and dataset names are unique")
			
	
	# make dataframe
	toDo = pd.DataFrame(rows, columns=['dataset', 'config_dataset', 'sample', 'host', 'host_fasta', 'virus', 'virus_fasta', 'merge', 'dedup', 'unique', 'outdir', 'bwa_mem_params', 'R1_file', 'R2_file', 'bam_file', 'adapter_1', 'adapter_2', 'postargs'])
	
	# do checks on dataframe
	check_dataset_sample_unique(toDo)
	
	ref_names = make_reference_dict(toDo)
	check_fastas_unique(toDo, ref_names)
	
	return toDo


def check_dataset_sample_unique(toDo):
	# check that every combination of 'dataset' and 'sample' is unique
	if len(set(toDo.loc[:,'unique'])) != len(toDo.loc[:,'unique']):
		raise InputError("Every combination of 'dataset' and 'sample' must be unique! Check inputs and try again")

def check_fastas_unique(toDo, ref_names):
	# check that each host/virus name refers to only one fasta
	for i, virus_name in enumerate(toDo.loc[:,'virus']):
		if toDo.loc[i,'virus_fasta'] != ref_names[virus_name]:
			raise InputError(f"Virus {virus_name} is used as a name in multiple datasets for different fasta files.  Please modify your configfile so that each unique virus/host name only refers to one fasta file")
	for i, host_name in enumerate(toDo.loc[:,'host']):
		if toDo.loc[i,'host_fasta'] != ref_names[host_name]:
			raise InputError(f"Host {host_name} is used as a name in multiple datasets for different fasta files.  Please modify your configfile so that each unique virus/host name only refers to one fasta file")
	
	
def make_reference_dict(toDo):
	# construct dictionary with reference names as keys and reference fastas as values
	
	ref_names = {}
	for index, row in toDo.iterrows():
		ref_names[row['host']] = row['host_fasta']
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
		
