#### pipeline for detecting viral integrations in NGS data ####

# the pipeline detecets viral integrations into a host genome in short reads by 
# identifying chimeric reads, discordant read pairs, and short insertions of virus 
# flanked by host sequence on both sides

# note that only paired-end reads are supported
#
# A yaml file specifies which references are to be used with which datasets
# the pipline lives in the intvi_pipeline folder which should be cloned
# from the git repo

# the steps in the pipeline are as follows:
# 1. preprocessing
#   - convert bam files to fastq, if necessary
#   - either merge R1 and R2 using seqprep, or don't do merging 
#     and just concatenate R1 and R2, as specified in the dsets.yaml config file
# 2. alignments
#   - align all the reads to the virus
#   - remove duplicates
#   - extract only the alinged reads to a fastq
#   - align these reads to the host
# 3. perl scripts
#   - run perl scripts with alignment sam  files as inputs to detect integrations
# 4. postprocessing
#   - apply various types of postprocessing: dedup, filter, mask, annotate
#   - generate ouput xlsx files for each dataset
#   - write output files for visualzation in UCSC browser
#  
# the primary output of the pipeline is the xlsx files for each datasets which
# contain the details of the detected integrations.


#### python modules ####

from glob import glob
from os import path, getcwd
import pandas as pd


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

#### CONFIG FILE ####
# supply on command line

bwa_mem_params="-A 1 -B 2 -O 6,6 -E 1,1 -L 0,0 -T 10 -h 200"


#### check file extensions ####

# supported file extensions for fastq files are .fastq and .fq, and compression with .gz and .bz2
# supported file extensions for bam/sam files are .bam, .sam

for dataset in config:

	if "R1_suffix" in config[dataset]:
		# check R1_suffix for compression and fastq file extension
		first_extension = path.splitext(config[dataset]["R1_suffix"])[1]
		remaining_suffix = path.splitext(config[dataset]["R1_suffix"])[0]
		# if uncompressed
		if first_extension == ".fq" or first_extension == ".fastq":
			# check read2 extension is the same as read1 extension
			if "R2_suffix" not in config[dataset]:
				raise InputError("If R1_suffix is specified (for fastq input), R2_suffix must also be specified")
			else:
				second_extension = path.splitext(config[dataset]["R2_suffix"])[1]
				if first_extension != second_extension:
					raise InputError("R1 and R2 file extensions must match")
		#bz2 compression
		elif first_extension == ".bz2":
			if "R2_suffix" not in config[dataset]:
				raise InputError("If R1_suffix is specified (for fastq input), R2_suffix must also be specified")
			# check read2 extension is the same as read1 extension
			else:
				second_extension = path.splitext(config[dataset]["R2_suffix"])[1]
				if first_extension != second_extension:
					raise InputError("R1 and R2 file extensions must match")
		#gzip compression
		elif first_extension == ".gz":
			if "R2_suffix" not in config[dataset]:
				raise InputError("If R1_suffix is specified (for fastq input), R2_suffix must also be specified")
			# check read2 extension is the same as read1 extension
			else:
				second_extension = path.splitext(config[dataset]["R2_suffix"])[1]
				if first_extension != second_extension:
					raise InputError("R1 and R2 file extensions must match")
		# unrecgonised file extension
		else:
			raise InputError("Only uncompressed ('.fq', '.fastq'), bzip2 ('.bz2') or gzip ('.gz') fastq files are currently supported")
	elif "bam_suffix" in config[dataset]:
		extension = path.splitext(config[dataset]["bam_suffix"])[1]
		if extension != ".bam" and extension != ".sam":
			raise InputError("For aligned input, only '.bam' and '.sam' files are currently supported")
	else:
		raise InputError("Please specify either 'R1_suffix' and 'R2_suffix' or 'bam_suffix' in the config file")
		



#### get information for wildcards ####

# get the name of each sample in each dataset, and save information about 
# how to process it in a dataframe

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
		suffix = config[dataset]["R1_suffix"]
		config[dataset]['read_folder'] = path.normpath(config[dataset]['read_folder'])
		folder = config[dataset]['read_folder']	
		samples = [path.basename(f)[:-len(suffix)] for f in glob(f"{folder}/*{suffix}")]
		if len(samples) == 0:
			print(f"warning: no files found for dataset {dataset}")
	# get bam/sam files for input
	elif "bam_suffix" in config[dataset]:
		suffix = config[dataset]["bam_suffix"]
		folder = config[dataset]['read_folder']	
		config[dataset]["R1_suffix"] = "_1.fq.gz"
		config[dataset]["R2_suffix"] = "_2.fq.gz"
		samples = [path.basename(f)[:-len(suffix)] for f in glob(f"{folder}/*{suffix}")]
		if len(samples) == 0:
			print(f"warning: no files found for dataset {dataset}")
	else:
		raise InputError(f"please specify either bam_suffix or R1_suffix and R2_suffix for dataset {dataset}")
			
	# figure out if 'dedup' and 'merge' are true or false
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
	
	
	for sample in samples:
		rows.append((dataset, sample, config[dataset]["host_name"], config[dataset]["host_fasta"], config[dataset]["virus_name"], config[dataset]["virus_fasta"], merge, dedup, f"{dataset}+++{sample}", config[dataset]['out_dir']))

toDo = pd.DataFrame(rows, columns=['dataset', 'sample', 'host', 'host_fasta', 'virus', 'virus_fasta', 'merge', 'dedup', 'unique', 'outdir'])

# check that every combination of 'dataset' and 'sample' is unique
if len(set(toDo.loc[:,'unique'])) != len(toDo.loc[:,'unique']):
	raise InputError("Every combination of 'dataset' and 'sample' must be unique! Check inputs and try again")

# construct dictionary with reference names as keys and reference fastas as values
host_names = { config[dataset]['host_name']:config[dataset]['host_fasta'] for dataset in config }
virus_names = { config[dataset]['virus_name']:config[dataset]['virus_fasta'] for dataset in config }
ref_names = {**host_names, **virus_names}

# check that each host/virus name refers to only one fasta
for i, virus_name in enumerate(toDo.loc[:,'virus']):
	if toDo.loc[i,'virus_fasta'] != ref_names[virus_name]:
		raise InputError(f"Virus {virus_name} is used as a name in multiple datasets for different fasta files.  Please modify your configfile so that each unique virus/host name only refers to one fasta file")
for i, host_name in enumerate(toDo.loc[:,'host']):
	if toDo.loc[i,'host_fasta'] != ref_names[host_name]:
		raise InputError(f"Host {host_name} is used as a name in multiple datasets for different fasta files.  Please modify your configfile so that each unique virus/host name only refers to one fasta file")

# construct arguments for postprocess.R script for each dataset
POSTARGS = {}
TOSORT = []
SORTED = []
for dataset in config:
	POSTARGS[dataset] = []
	if "post" in config[dataset]:
		for element in config[dataset]["post"]:
			# need to check if this element is a string or a dict
			if type(element) is str:
				# look for keys to be 'filter' or 'dedup'
				if element == "filter":
					POSTARGS[dataset].append("filter")
				elif element == "dedup":
					POSTARGS[dataset].append("dedup")
			# postprocessing types with files specified will be in ordered dictionaries
			elif type(element) is OrderedDict:
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


#### global wildcard constraints ####

wildcard_constraints:
	virus="|".join(set(toDo.loc[:,'virus'])),
	samp="|".join(set(toDo.loc[:,'sample'])),
	dset="|".join(set(toDo.loc[:,'dataset'])),
	host="|".join(set(toDo.loc[:,'host'])),
	align_type="bwaPaired|bwaSingle"


#### local rules ####
localrules: all, combine, check_bam_input_is_paired

#### target files ####
rule all:
	input: 
		#"{outpath}/summary/count_mapped.txt",
		expand("{outpath}/summary/{dset}.xlsx", zip, outpath = toDo.loc[:, 'outdir'], dset=toDo.loc[:,'dataset']),
		expand("{outpath}/summary/ucsc_bed/{dset}.post.bed", zip, outpath = toDo.loc[:, 'outdir'], dset=toDo.loc[:,'dataset'])

#### read preprocessing ####

rule check_bam_input_is_paired:
	input: lambda wildcards: path.normpath(f"{config[wildcards.dset]['read_folder']}/{wildcards.samp}{config[wildcards.dset]['bam_suffix']}")
	output:
		ok = temp("{outpath}/{dset}/.reads/{samp}.tmp"),
	conda:
		"envs/bwa.yml"
	shell:
		"""
		FWD=$(samtools view -c -f 0x40 {input})
		if [[ "$FWD" == "0" ]]; then
			echo "Data must be paired-end"
			exit 1
		fi
		REV=$(samtools view -c -f 0x80 {input})
		if [[ "$FWD" != "$REV" ]]; then
			echo "Number of forward reads must match number of paired reads"
			exit 1
		fi
		touch {output.ok}
		"""	

rule bam_to_fastq:
	input:
		bam = lambda wildcards: path.normpath(f"{config[wildcards.dset]['read_folder']}/{wildcards.samp}{config[wildcards.dset]['bam_suffix']}"),
		ok = lambda wildcards: f"{wildcards.outpath}/{wildcards.dset}/.reads/{wildcards.samp}.tmp"
	output:
		sorted_bam = temp("{outpath}/{dset}/.reads/{samp}.sorted.bam"),
		r1 = temp("{outpath}/{dset}/.reads/{samp}_1.fq.gz"),
		r2 = temp("{outpath}/{dset}/.reads/{samp}_2.fq.gz"),
	conda:
		"envs/picard.yml"
	shell:
		"""
		picard SortSam I={input.bam} O={output.sorted_bam} SORT_ORDER=queryname
		picard SamToFastq I={output.sorted_bam} FASTQ={output.r1} SECOND_END_FASTQ={output.r2}
		"""


# input functions for if had a bam or fastq as input
def get_for_seqprep(wildcards, read_type):
	row_idx = list(toDo.loc[:,'unique']).index(f"{wildcards.dset}+++{wildcards.samp}")
	# get true/True values for dedup
	if 'bam_suffix' in config[wildcards.dset]:
		if read_type == "1":
			return "{outpath}/{dset}/.reads/{samp}_1.fq.gz"
		else:
			return "{outpath}/{dset}/.reads/{samp}_2.fq.gz"
	else:
		if read_type == "1":
			return path.normpath(f"{config[wildcards.dset]['read_folder']}/{wildcards.samp}{config[wildcards.dset]['R1_suffix']}")
		else:
			return path.normpath(f"{config[wildcards.dset]['read_folder']}/{wildcards.samp}{config[wildcards.dset]['R2_suffix']}")


rule seqPrep:
# if we're doing it
	input:
		r1 = lambda wildcards: get_for_seqprep(wildcards, "1"),
		r2 = lambda wildcards: get_for_seqprep(wildcards, "2")
	output:
		merged = temp("{outpath}/{dset}/.merged_reads/{samp}.SeqPrep_merged.fastq.gz"),
		proc_r1 = temp("{outpath}/{dset}/.merged_reads/{samp}.1.fastq.gz"),
		proc_r2 = temp("{outpath}/{dset}/.merged_reads/{samp}.2.fastq.gz"),
		all = temp("{outpath}/{dset}/.merged_reads/{samp}.all.fastq.gz")
	conda:	
		"envs/seqprep.yml"
	params:
		A = lambda wildcards: config[wildcards.dset]["read1-adapt"],
		B = lambda wildcards: config[wildcards.dset]["read2-adapt"]
	shell:
		"""
		SeqPrep -A {params.A} -B {params.B} -f {input.r1} -r {input.r2} -1 {output.proc_r1} -2 {output.proc_r2} -s {output.merged}
		cat {output.proc_r1} {output.proc_r2} {output.merged} > {output.all}
		"""

# helper function to figure out which cat to use for rule combine (depending on compression)
def get_compression(input):
	ext = path.splitext(input.r1)[1]
	if ext == ".bz2":
		return "bzcat"
	elif ext == ".gz":
		return "zcat"
	else:
		return "cat"


rule combine:
# for if we don't want to do seqPrep
	input:
		r1 = lambda wildcards: get_for_seqprep(wildcards, "1"),
		r2 = lambda wildcards: get_for_seqprep(wildcards, "2")
	output:
		proc_r1 = temp("{outpath}/{dset}/.combined_reads/{samp}.1.fastq.gz"),
		proc_r2 = temp("{outpath}/{dset}/.combined_reads/{samp}.2.fastq.gz"),
		all = temp("{outpath}/{dset}/.combined_reads/{samp}.all.fastq.gz")
	params:
		test = lambda wildcards, input: print(get_compression(input)),
		cat = lambda wildcards, input: get_compression(input),
		
	shell:
		"""
		{params.cat} {input.r1} | gzip > {output.proc_r1}
		{params.cat} {input.r2} | gzip > {output.proc_r2}
		{params.cat} {input.r1} {input.r2} | gzip > {output.all}
		"""
		
#functions for if we did seqPrep or not
def get_for_align(wildcards, read_type):
	row_idx = list(toDo.loc[:,'unique']).index(f"{wildcards.dset}+++{wildcards.samp}")
	merge_check = bool(toDo.loc[toDo.index[row_idx], 'merge'])
	if merge_check is True:
		folder = ".merged_reads"
	else:
		folder = ".combined_reads"
	return f"{wildcards.outpath}/{wildcards.dset}/{folder}/{wildcards.samp}.{read_type}.fastq.gz"


#### alignments ####

rule index:
	input:
		lambda wildcards: ref_names[wildcards.genome]
	output:
		expand("{outpath}/.references/{genome}/{genome}.{ext}", ext=["ann", "amb", "bwt", "pac", "sa"], allow_missing=True)
	conda: 
		"envs/bwa.yml"
	params:
		prefix = lambda wildcards, output: path.splitext(output[0])[0]
	shell:
		"bwa index -p {params.prefix} {input}"

	
rule align_bwa_virus_single:
	input:
		idx = expand("{outpath}/.references/{virus}/{virus}.{ext}", ext=["ann", "amb", "bwt", "pac", "sa"], allow_missing=True),
		all = lambda wildcards: get_for_align(wildcards, "all"),
	output:
		sam = temp("{outpath}/{dset}/.virus_aligned/{samp}.{virus}.bwaSingle.sam")
	params:
		index = lambda wildcards, input: path.splitext(input.idx[0])[0],
		mapping = bwa_mem_params
	conda:
		"envs/bwa.yml"
	threads: 5
	shell:
		"""
		bwa mem -t {threads} {params.mapping} -o {output.sam} {params.index} {input.all}
		"""

rule align_bwa_virus_paired:
	input:
		idx = expand("{outpath}/.references/{virus}/{virus}.{ext}", ext=["ann", "amb", "bwt", "pac", "sa"], allow_missing=True),
		r1 = lambda wildcards: get_for_align(wildcards, "1"),
		r2 = lambda wildcards: get_for_align(wildcards, "2"),
	output:
		sam = temp("{outpath}/{dset}/.virus_aligned/{samp}.{virus}.bwaPaired.sam")
	params:
		index = lambda wildcards, input: path.splitext(input.idx[0])[0],
		mapping = bwa_mem_params
	conda:
		"envs/bwa.yml"
	threads: 5
	shell:
		"""
		bwa mem -t {threads} {params.mapping} -o {output.sam} {params.index} {input.r1} {input.r2} 
		"""
		

def get_sam(wildcards, readType, genome):
	row_idx = list(toDo.loc[:,'unique']).index(f"{wildcards.dset}+++{wildcards.samp}")
	dedup_check = bool(toDo.loc[toDo.index[row_idx], 'dedup'])
	# if we want host alignment
	if genome == "virus":
		# if we're doing deduplication
		if dedup_check is True:
			# if we want single reads
			if readType == "single":
				return path.splitext(rules.align_bwa_virus_single.output[0])[0] + ".rmdup.sam"
			else:
				return path.splitext(rules.align_bwa_virus_paired.output[0])[0] + ".rmdup.sam"
		# if we're not doing deduplication
		else:
			if readType == "single":
				return path.splitext(rules.align_bwa_virus_single.output[0])[0] + ".dups.sam"
			else:
				return path.splitext(rules.align_bwa_virus_paired.output[0])[0] + ".dups.sam"
	# if we want the host alignment
	else:
		# if we're doing deduplication
		if dedup_check is True:
			# if we want single reads
			if readType == "single":
				return path.splitext(rules.align_bwa_host_single.output[0])[0] + ".rmdup.sam"
			else:
				return path.splitext(rules.align_bwa_host_paired.output[0])[0] + ".rmdup.sam"
		# if we're not doing deduplication
		else:
			if readType == "single":
				return path.splitext(rules.align_bwa_host_single.output[0])[0] + ".dups.sam"
			else:
				return path.splitext(rules.align_bwa_host_paired.output[0])[0] + ".dups.sam"	
				
		
		
rule extract_vAligned_single:
	input:
		aligned = lambda wildcards: get_sam(wildcards, "single", "virus"),
	output:
		sam = temp("{outpath}/{dset}/.virus_aligned/{samp}.{virus}.bwaSingle.mapped.sam"),
	conda:
		"envs/bwa.yml"
	shell:
		"""
		samtools view -h -F 0x4 -F 0x800 -o - {input.aligned} | samtools sort - -n -o {output.sam}
		"""

rule extract_vAligned_paired:
	input:
		aligned = lambda wildcards: get_sam(wildcards, "paired", "virus")
	output:
		pvBam_readMap_mateUnmap = temp("{outpath}/{dset}/.virus_aligned/{samp}.{virus}.bwaPaired.L.bam"),
		pvBam_readUnmap_mateMap = temp("{outpath}/{dset}/.virus_aligned/{samp}.{virus}.bwaPaired.R.bam"),
		pvBam_bothMapped = temp("{outpath}/{dset}/.virus_aligned/{samp}.{virus}.bwaPaired.B.bam"),
		sam = temp("{outpath}/{dset}/.virus_aligned/{samp}.{virus}.bwaPaired.mapped.sam"),
	conda:
		"envs/bwa.yml"
	shell:
		"""
		samtools view -hb -F 0x4 -f 0x8 -F 0x800 -o {output.pvBam_readMap_mateUnmap} {input.aligned}
		samtools view -hb -f 0x4 -F 0x8 -F 0x800 -o {output.pvBam_readUnmap_mateMap} {input.aligned}
		samtools view -hb -F 0x4 -F 0x8 -F 0x800 -o {output.pvBam_bothMapped} {input.aligned}
		samtools merge - {output.pvBam_readMap_mateUnmap} {output.pvBam_bothMapped} {output.pvBam_readUnmap_mateMap} | samtools sort - -n -o {output.sam}
		"""
		

rule extract_to_fastq_single:
	input:
		sam = rules.extract_vAligned_single.output.sam
	output:
		fastq = temp("{outpath}/{dset}/.virus_aligned/{samp}.bwaSingle.mappedTo{virus}.fastq.gz")
	conda:
		"envs/picard.yml"
	shell:
		"""
		picard SamToFastq I={input.sam} FASTQ={output.fastq}
		"""
		
rule extract_to_fastq_paired:
	input:
		sam = rules.extract_vAligned_paired.output.sam
	output:
		fastq1 = temp("{outpath}/{dset}/.virus_aligned/{samp}.bwaSingle.mappedTo{virus}.1.fastq.gz"),
		fastq2 = temp("{outpath}/{dset}/.virus_aligned/{samp}.bwaSingle.mappedTo{virus}.2.fastq.gz")
	conda:
		"envs/picard.yml"
	shell:
		"""
		picard SamToFastq I={input.sam} FASTQ={output.fastq1} SECOND_END_FASTQ={output.fastq2}
		"""

rule align_bwa_host_single:
	input:	
		idx = expand("{outpath}/.references/{host}/{host}.{ext}", ext=["ann", "amb", "bwt", "pac", "sa"], allow_missing=True),
		all = rules.extract_to_fastq_single.output.fastq
	output:
		sam = temp("{outpath}/{dset}/.host_aligned/{samp}.{host}.readsFrom{virus}.bwaSingle.sam"),
	conda: 
		"envs/bwa.yml"
	params:
		index = lambda wildcards, input: path.splitext(input.idx[0])[0],
		mapping = bwa_mem_params
	threads: 4
	shell:		
		"""
		bwa mem -t {threads} {params.mapping} -o {output.sam} {params.index} {input.all}
		"""
		
rule align_bwa_host_paired:
	input:	
		idx = expand("{outpath}/.references/{host}/{host}.{ext}", ext=["ann", "amb", "bwt", "pac", "sa"], allow_missing=True),
		r1 = rules.extract_to_fastq_paired.output[0],
		r2 = rules.extract_to_fastq_paired.output[1]
	output:
		sam = temp("{outpath}/{dset}/.host_aligned/{samp}.{host}.readsFrom{virus}.bwaPaired.sam"),
	conda: 
		"envs/bwa.yml"
	params:
		index = lambda wildcards, input: path.splitext(input.idx[0])[0],
		mapping = bwa_mem_params
	threads: 4
	shell:		
		"""
		bwa mem -t {threads} {params.mapping} -o {output.sam} {params.index} {input.r1} {input.r2}
		"""

#### sam file manipulations ####

rule convert_to_bam:
	input:
		"{outpath}/{dset}/{folder}/{alignment}.sam"
	output:
		bam = "{outpath}/{dset}/{folder}/{alignment}.bam",
		bai = "{outpath}/{dset}/{folder}/{alignment}.bam.bai"
	params:
		tmp_prefix = lambda wildcards, input: path.splitext(input[0])[0]
	wildcard_constraints:
		folder = "\.host_aligned|\.virus_aligned"
	conda: 
		"envs/bwa.yml"	
	shell:
		"""
		rm -f {params.tmp_prefix}*tmp*
		samtools view -bhS {input} | samtools sort - -o {output.bam}
		samtools index {output.bam}
		"""
		
rule markdup:
	input:
		sam = "{outpath}/{dset}/{folder}/{alignment}.sam"
	output:
		fixmate = temp("{outpath}/{dset}/{folder}/{alignment}.fixmate.bam"),
		markdup = "{outpath}/{dset}/{folder}/{alignment}.dups.sam",
		metrics = temp("{outpath}/{dset}/{folder}/{alignment}.dups.txt")
	wildcard_constraints:
		folder = "\.host_aligned|\.virus_aligned"
	conda: 
		"envs/picard.yml"	
	shell:
		"""
		picard FixMateInformation I={input.sam} O={output.fixmate} ADD_MATE_CIGAR=true SORT_ORDER=queryname
		picard MarkDuplicates I={output.fixmate} O={output.markdup} METRICS_FILE={output.metrics}
		"""
		
rule rmdup:
	input:
		sam = "{outpath}/{dset}/{folder}/{alignment}.dups.sam"	
	output:
		sam = temp("{outpath}/{dset}/{folder}/{alignment}.rmdup.sam")	
	wildcard_constraints:
		folder = "\.host_aligned|\.virus_aligned"
	conda: 
		"envs/bwa.yml"	
	shell:
		"""
		samtools view -h -F 1024 {input.sam} > {output.sam}
		"""
	

#### perl scripts ####

rule run_soft:
	input:
		host = lambda wildcards: get_sam(wildcards, "single", "host"),
		virus = rules.extract_vAligned_single.output[0],
	output:
		soft = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.soft.txt"),
	shell:
		"""
		perl -I. ./softClip.pl --viral {input.virus} --human {input.host} --output {output.soft} --tol 3
		"""
		
rule run_short:
	input:
		host = lambda wildcards: get_sam(wildcards, "single", "host"),
		virus = rules.extract_vAligned_single.output[0],
	output:
		short = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.short.txt"),
	shell:
		"""
		perl -I. ./short.pl --viral {input.virus} --human {input.host} --output {output.short} --tol 3
		"""
		
rule run_discordant:
	input:
		host = lambda wildcards: get_sam(wildcards, "paired", "host"),
		virus = rules.extract_vAligned_paired.output[3],
	output:
		discord = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.discordant.txt"),
	shell:
		"""
		perl -I. ./discordant.pl --viral {input.virus} --human {input.host} --output {output.discord} --tol 3
		"""

rule combine_ints:
	input:
		soft = rules.run_soft.output,
		short = rules.run_short.output,
		discordant = rules.run_discordant.output
	output:
		temp =  temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.txt.tmp"),
		all = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.txt"
	shell:
		"""
		awk 'FNR>1 || NR==1' {input} > {output.all}
		sort -n -k1,1 -k2,2n {output.all} > {output.temp}
		cp {output.temp} {output.all}
		"""

	
#### postprocessing ####

rule sortbed:
	input:
		TOSORT
	output:
		SORTED
	run:
		for i in range(len(TOSORT)):
			print(f"sorting file {TOSORT[i]} into file {SORTED[i]}, file extension {path.splitext(TOSORT[i])[1]}")
			
			# if file is a bed file
			if (path.splitext(TOSORT[i])[1] == ".bed"):
				body = f"sort -k1,1 -k2,2n {TOSORT[i]} > {SORTED[i]}"
				shell(body)
				print(body)
			# if file is a gtf file
			elif (path.splitext(TOSORT[i])[1] == ".gtf"):
				
				body = f"awk '{{{{ if ($0 !~ /^#/) {{{{ print $0 }}}} }}}}' {TOSORT[i]} | sort -k1,1 -k4,4n > {SORTED[i]}"
				shell(body)
				print(body)
			else:
				raise ValueError("only gtf files and bed files are supported")


#this rule (post) sometimes causes issues with conda. The error is usually something to do with ldpaths:

#Activating conda environment: /scratch1/sco305/intvi_cmri/intvi_pipeline/.snakemake/conda/586e76e5
#/scratch1/sco305/intvi_cmri/intvi_pipeline/.snakemake/conda/586e76e5/lib/R/bin/R: line 238: /scratch1/sco305/intvi_cmri/#intvi_pipeline/.snakemake/conda/586e76e5/lib/R/etc/ldoutpaths: No such file or directory

# however, sometimes this rule runs just fine.

# this issue is descrived here: https://github.com/conda-forge/r-base-feedstock/issues/67
# however, it doesn't appear to have been resolved.  temporarily get around this by re-attempting jobs
# that failed using command line option --restart-times, but will need to come up with a better solution for this

rule post:
	input:
		"{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.txt",
		rules.sortbed.output
	output:
		"{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.post.txt"
	conda:
		"envs/rscripts.yml"
	params:
		lambda wildcards, input: f"{input[0]} {POSTARGS[wildcards.dset]}"
	shell:
		"""
		Rscript post/postprocess.R {params}
		"""

	
rule summarise:
	input:
		lambda wildcards: [f"{wildcards.outpath}/{wildcards.dset}/ints/{samp}.{host}.{virus}.integrations.post.txt" for samp, host, virus
			in zip(toDo.loc[toDo['dataset'] == wildcards.dset,'sample'], 
					toDo.loc[toDo['dataset'] == wildcards.dset,'host'], 
					toDo.loc[toDo['dataset'] == wildcards.dset,'virus'])]
	output:
		"{outpath}/summary/{dset}.xlsx",
		"{outpath}/summary/{dset}_annotated.xlsx"
	conda:
		"envs/rscripts.yml"
	shell:
		"Rscript summarise_ints.R {input}"

rule ucsc_bed:
	input:
		expand("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.post.txt", 
					zip, 
					outpath = toDo.loc[:,'outdir'],
					dset = toDo.loc[:,'dataset'], 
					samp = toDo.loc[:,'sample'], 
					host = toDo.loc[:,'host'], 
					virus = toDo.loc[:,'virus'], 
					allow_missing=True),
	output:
		expand("{outpath}/summary/ucsc_bed/{dset}.post.bed", zip, outpath = toDo.loc[:, 'outdir'], dset = toDo.loc[:,'dataset']),
	conda:
		"envs/rscripts.yml"
	shell:
		"""
		Rscript writeBed.R {input}
		bash -e format_ucsc.sh
		"""
