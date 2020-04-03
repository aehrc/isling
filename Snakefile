#### pipeline for detecting viral integrations in NGS data ####

# the pipeline detecets viral integrations into a host genome in short reads by 
# identifying chimeric reads, discordant read pairs, and short insertions of virus 
# flanked by host sequence on both sides

# note that only paired-end reads are supported

# the pipeline relies on a particular directory structure:
 
# parent_folder/
# ├── data
# │   ├── metadata
# │   │   └── dataset1
# │   ├── reads
# │   │   └── dataset1
# │   │       ├── sample1_1.fastq.gz
# │   │       ├── sample2_2.fastq.gz
# │   │       ├── sample2_1.fastq.gz
# │   │       └── sample2_2.fastq.gz
# │   └── references
# │       ├── host.fa
# │       └── virus.fa
# ├── intvi_pipeline 
# │   └── all pipeline files from git repo
# └── proj-spcfc
#     └── dsets.yaml
#
# the folder data contains all the data necessary to run the pipeline
# within the data folder, reads should be organised into folders representing
# datasets, with gzipped paired fastq files for each sample named as above
# references should be placed in data/references.  A yaml file in proj-spcfc
# specifies which references are to be used with which datasets
# the pipline lives in the intvi_pipeline folder which should be cloned
# from the git repo

# the steps in the pipeline are as follows:
# 1. preprocessing
#   - rename all files for consistency (files may either be named "<sample>_1.fastq.gz" or "<sample>_R1.fastq.gz")
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
from os import path
import pandas as pd

#### CONFIG FILE ####
# supply on command line

references_path = "../data/references/"

#make four lists, with each entry in each list corresponding to one desired output file

DATASETS = []
SAMPLES = []
HOSTS = []
VIRUSES = []
#also make dictionary MERGE, with key dataset_sample and value True or False
MERGE = {}

rows = []

# we need to know if input files are bam or fastq.  check for filename extensions
# currently only support fastq (.fq.gz, .fq, .fastq.gz, .fastq)
# and bam (.bam, .sam)

# get files for each dataset, and save all paramters in a pandas dataframe
for dataset in config:
	# for fastq files
	if "R1_suffix" in config[dataset]:
		suffix = config[dataset]["R1_suffix"]	
		samples = [path.basename(f)[:-len(suffix)] for f in glob(f"../data/reads/{dataset}/*{suffix}")]
		if len(samples) == 0:
			print(f"warning: no files found for dataset {dataset}")
		config[dataset]['bam_suffix'] = ""
	# for bam/sam files
	elif "bam_suffix" in config[dataset]:
		suffix = config[dataset]["bam_suffix"]
		config[dataset]["R1_suffix"] = "_1.fq"
		config[dataset]["R2_suffix"] = "_2.fq"
		samples = [path.basename(f)[:-len(suffix)] for f in glob(f"../data/reads/{dataset}/*{suffix}")]
		if len(samples) == 0:
			print(f"warning: no files found for dataset {dataset}")
	else:
		raise InputError(f"please specify either bam_suffix or R1_suffix and R2_suffix for dataset {dataset}")
			
	for sample in samples:
		rows.append((dataset, sample, config[dataset]["host"], config[dataset]["virus"], config[dataset]["merge"], config[dataset]["dedup"], f"{dataset}+++{sample}"))

toDo = pd.DataFrame(rows, columns=['dataset', 'sample', 'host', 'virus', 'merge', 'dedup', 'unique'])

# check that every combination of 'dataset', 'sample', 'host' and 'virus' is unique
if len(set(toDo.loc[:,'unique'])) != len(toDo.loc[:,'unique']):
	raise InputError("Every combination of 'dataset' and 'sample' must be unique! Check inputs and try again")

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


#### global wildcard constraints ####

wildcard_constraints:
	virus="|".join(set(toDo.loc[:,'virus'])),
	sample="|".join(set(toDo.loc[:,'sample'])),
	dset="|".join(set(toDo.loc[:,'dataset'])),
	host="|".join(set(toDo.loc[:,'host'])),
	align_type="bwaShort|bwaPaired|bwaSoft"


#### local rules ####
localrules: all, combine, 

#### target files ####
rule all:
	input: 
		#"../out/summary/count_mapped.txt",
		expand("../out/summary/{dset}.xlsx", dset=set(toDo.loc[:,'dataset'])),
		#expand("../out/summary/ucsc_bed/{dset}.post.bed", dset=set(toDo.loc[:,'dataset'])),
		#expand("../out/{dset}/ints/{samp}.{host}.{virus}.integrations.post.txt", 
		#	zip, 
		#	dset=toDo.loc[:,'dataset'],
		#	samp=toDo.loc[:,'sample'], 
		#	host=toDo.loc[:,'host'], 
		#	virus=toDo.loc[:,'virus'])

#### read preprocessing ####

rule bamtobed:
	input: lambda wildcards: f"../data/reads/{wildcards.dset}/{wildcards.samp}{config[wildcards.dset]['bam_suffix']}"
	output: 
		r1 = temp("../data/reads/{dset}/{samp}_1.fq"),
		r2 = temp("../data/reads/{dset}/{samp}_2.fq"),
	params:
		sortedBam = lambda wildcards: f"../out/reads/{wildcards.dset}/{wildcards.samp}.sorted{config[wildcards.dset]['bam_suffix']}",
		unzip_r1 = lambda wildcards: f"../data/reads/{wildcards.dset}/{wildcards.samp}_1.fastq",
		unzip_r2 = lambda wildcards: f"../data/reads/{wildcards.dset}/{wildcards.samp}_2.fastq",
		tmpfiles = lambda wildcards, output: f"{path.dirname(output.r1)}/*.tmp"
	conda:
		"envs/bwa.yml"
	shell:
		"""
		rm {params.tmpfiles}
		dirname {params.sortedBam} | xargs mkdir -p
		samtools sort -n -o {params.sortedBam} {input}
		bedtools bamtofastq -i {params.sortedBam} -fq {params.unzip_r1} -fq2 {params.unzip_r2}
		rm {params.sortedBam}
		"""

rule dedup:
# remove exact duplicates (subs=0)
	input:
		r1 = lambda wildcards: f"../data/reads/{wildcards.dset}/{wildcards.samp}{config[wildcards.dset]['R1_suffix']}",
		r2 = lambda wildcards: f"../data/reads/{wildcards.dset}/{wildcards.samp}{config[wildcards.dset]['R2_suffix']}"
	output:
		r1 = temp("../out/{dset}/dedup_reads/{samp}.1.fastq.gz"),
		r2 = temp("../out/{dset}/dedup_reads/{samp}.2.fastq.gz"),
		all = temp("../out/{dset}/dedup_reads/{samp}.all.fastq.gz")
	conda:
		"envs/bbmap.yml"
	shell:
		"""
		clumpify.sh -Xmx15g in1="{input.r1}" in2="{input.r2}" out1="{output.r1}" out2="{output.r2}" dedupe subs=0
		cat {output.r1} {output.r2} > {output.all}
		"""

# input functions for if we want to do deduplication or not
def get_for_seqprep(wildcards, read_type):
	key = f"{wildcards.dset}+++{wildcards.samp}"
	# get row index for this key in toDo
	row_idx = list(toDo.loc[:,'unique']).index(key)
	# get true/True values for dedup
	dedup_check = toDo.iloc[row_idx, 5]
	if (dedup_check == "True") or (dedup_check == "true") or (dedup_check is True):
		return f"../out/{wildcards.dset}/dedup_reads/{wildcards.samp}.{read_type}.fastq.gz"
	else:
		if read_type == "1":
			return f"../data/reads/{wildcards.dset}/{wildcards.samp}{config[wildcards.dset]['R1_suffix']}"
		else:
			return f"../data/reads/{wildcards.dset}/{wildcards.samp}{config[wildcards.dset]['R2_suffix']}"


rule seqPrep:
# if we're doing it
	input:
		r1 = lambda wildcards: get_for_seqprep(wildcards, "1"),
		r2 = lambda wildcards: get_for_seqprep(wildcards, "2")
	output:
		merged = temp("../out/{dset}/merged_reads/{samp}.SeqPrep_merged.fastq.gz"),
		proc_r1 = temp("../out/{dset}/merged_reads/{samp}.1.fastq.gz"),
		proc_r2 = temp("../out/{dset}/merged_reads/{samp}.2.fastq.gz"),
		all = temp("../out/{dset}/merged_reads/{samp}.all.fastq.gz")
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

rule combine:
# for if we don't want to do seqPrep
	input:
		r1 = lambda wildcards: get_for_seqprep(wildcards, "1"),
		r2 = lambda wildcards: get_for_seqprep(wildcards, "2")
	output:
		proc_r1 = temp("../out/{dset}/combined_reads/{samp}.1.fastq.gz"),
		proc_r2 = temp("../out/{dset}/combined_reads/{samp}.2.fastq.gz"),
		all = temp("../out/{dset}/combined_reads/{samp}.all.fastq.gz")
	shell:
		"""
		cp {input.r1} {output.proc_r1}
		cp {input.r2} {output.proc_r2}
		cat {input.r1} {input.r2} > {output.all}
		"""
		
#functions for if we did seqPrep and deduplication or not
def get_for_align(wildcards, read_type):
	key = f"{wildcards.dset}+++{wildcards.samp}"
	# get row index for this key in toDo
	row_idx = list(toDo.loc[:,'unique']).index(key)
	# get true/True values for merge and dedup
	merge_check = toDo.iloc[row_idx, 4]
	dedup_check = toDo.iloc[row_idx, 5]
	if (merge_check == "True") or (merge_check == "true") or (merge_check is True):
		folder = "merged_reads"
	elif (dedup_check == "True") or (dedup_check == "true") or (dedup_check is True):
		folder = "dedup_reads"
	else:
		folder = "combined_reads"
	return f"../out/{wildcards.dset}/{folder}/{wildcards.samp}.{read_type}.fastq.gz"


#### alignments ####

rule index:
	input:
		"../data/references/{genome}.fa"
	output:
		expand("../out/.references/{genome}/{genome}.{ext}", ext=["ann", "amb", "bwt", "pac", "sa"], allow_missing=True)
	conda: 
		"envs/bwa.yml"
	params:
		prefix = lambda wildcards: references_path + wildcards.genome
	shell:
		"bwa index -p {params.prefix} {input}"

	
rule align_bwa_virus:
	input:
		idx = expand("../out/.references/{virus}/{virus}.{ext}", ext=["ann", "amb", "bwt", "pac", "sa"], allow_missing=True),
		all = lambda wildcards: get_for_align(wildcards, "all"),
	output:
		sam = temp("../out/{dset}/.virus_aligned/{samp}.{virus}.{align_type}.sam")
	params:
		index = lambda wildcards: references_path + wildcards.virus
	threads: 5
	shell:
		"""
		python ./alignReadsWithBWA.py --threads {threads} --index {params.index} --read1 {input.all} --output {output.sam} --threshold 10 --hflag 200 
		"""

		
rule extract_vAligned:
	input:
		vSing = rules.align_bwa_virus.output.sam,
	output:
		svSam = temp("../out/{dset}/.virus_aligned/{samp}.{virus}.{align_type}.mapped.sam"),
	conda:
		"envs/bwa.yml"
	shell:
		"""
		samtools view -h -F 0x4 -F 0x800 -o {output.svSam} {input.vSing} 
		"""


rule extract_to_fastq:
	input:
		sam = rules.extract_vAligned.output
	output:
		fastq = temp("../out/{dset}/.virus_aligned/{samp}.{align_type}.mappedTo{virus}.fastq.gz")
	conda:
		"envs/picard.yml"
	shell:
		"""
		picard SamToFastq I={input.sam} FASTQ={output.fastq}
		"""

rule align_bwa_host:
	input:	
		idx = expand("../out/.references/{host}/{host}.{ext}", ext=["ann", "amb", "bwt", "pac", "sa"], allow_missing=True),
		all = rules.extract_to_fastq.output
	output:
		sam = temp("../out/{dset}/.host_aligned/{samp}.{host}.readsFrom{virus}.{align_type}.sam"),
	conda: 
		"envs/bwa.yml"
	params:
		index = lambda wildcards: references_path + wildcards.host
	threads: 5
	shell:		
		"""
		python ./alignReadsWithBWA.py --threads {threads} --index {params.index} --read1 {input.all} --output {output.sam} --threshold 10 --hflag 200 --insert_extend 0
		"""

rule convert:
	input:
		"../out/{dset}/{host_virus}/{name}.sam"
	output:
		bam = "../out/{dset}/{host_virus}/{name}.bam",
		bai = "../out/{dset}/{host_virus}/{name}.bam.bai"
	conda: 
		"envs/bwa.yml"	
	shell:
		"""
		samtools view -bhS {input} | samtools sort - -o {output.bam}
		samtools index {output.bam}
		"""

#### perl scripts ####

rule run_soft:
	input:
		host ="../out/{dset}/.host_aligned/{samp}.{host}.readsFrom{virus}.bwaSoft.sam",
		virus = "../out/{dset}/.virus_aligned/{samp}.{virus}.bwaSoft.sam",
	output:
		soft = temp("../out/{dset}/ints/{samp}.{host}.{virus}.soft.txt"),
	shell:
		"""
		perl -I. ./softClip.pl --viral {input.virus} --human {input.host} --output {output.soft} --tol 3
		"""

rule combine_ints:
	input:
		soft = "../out/{dset}/ints/{samp}.{host}.{virus}.soft.txt"
	output:
		temp =  "../out/{dset}/ints/{samp}.{host}.{virus}.integrations.txt.tmp",
		all = "../out/{dset}/ints/{samp}.{host}.{virus}.integrations.txt"
	shell:
		"""
		sed -e '2,${{/Chr/d' -e '}}' {input} > {output.all}
		sort -n -k1,1 -k2,2n {output.all} > {output.temp}
		cp {output.temp} {output.all}
		"""

	
#### postprocessing ####

rule count_mapped:
	input:
		expand(expand("../out/{dset}/.virus_aligned/{samp}.{virus}.{align_type}.bam", 
					zip, 
					dset = toDo.loc[:,'dataset'], 
					samp = toDo.loc[:,'sample'], 
					virus = toDo.loc[:,'virus'],
					allow_missing=True),
						align_type=['bwaSoft']),
		#expand(expand("../out/{dset}/.host_aligned/{samp}.{host}.readsFrom{virus}.{align_type}.bam", 
		#			zip, 
		#			dset = toDo.loc[:,'dataset'], 
		#			samp = toDo.loc[:,'sample'], 
		#			virus = toDo.loc[:,'virus'],
		#			host = toDo.loc[:,'host'],
		#			allow_missing=True),
		#				align_type=['bwaSoft']),
						
	output:
		"../out/summary/count_mapped.txt"
	conda: 
		"envs/bwa.yml"
	shell:
		"./count_mapped.sh"

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
#/scratch1/sco305/intvi_cmri/intvi_pipeline/.snakemake/conda/586e76e5/lib/R/bin/R: line 238: /scratch1/sco305/intvi_cmri/#intvi_pipeline/.snakemake/conda/586e76e5/lib/R/etc/ldpaths: No such file or directory

# however, sometimes this rule runs just fine.

# this issue is descrived here: https://github.com/conda-forge/r-base-feedstock/issues/67
# however, it doesn't appear to have been resolved.  temporarily get around this by re-attempting jobs
# that failed using command line option --restart-times, but will need to come up with a better solution for this

rule post:
	input:
		"../out/{dset}/ints/{samp}.{host}.{virus}.integrations.txt",
		rules.sortbed.output
	output:
		"../out/{dset}/ints/{samp}.{host}.{virus}.integrations.post.txt"
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
		expand("../out/{dset}/ints/{samp}.{host}.{virus}.integrations.post.txt", 
			zip, 
			samp=toDo.loc[:,'sample'], 
			host=toDo.loc[:,'host'], 
			virus=toDo.loc[:,'virus'], 
			allow_missing=True)

	output:
		"../out/summary/{dset}.xlsx"
	conda:
		"envs/rscripts.yml"
	shell:
		"Rscript summarise_ints.R {input}"

rule ucsc_bed:
	input:
		expand("../out/{dset}/ints/{samp}.{host}.{virus}.integrations.post.txt", 
					zip, 
					dset = toDo.loc[:,'dataset'], 
					samp = toDo.loc[:,'sample'], 
					host = toDo.loc[:,'host'], 
					virus = toDo.loc[:,'virus']),
	output:
		expand("../out/summary/ucsc_bed/{dset}.post.bed", dset = set(toDo.loc[:,'dataset'])),
	conda:
		"envs/rscripts.yml"
	shell:
		"""
		Rscript writeBed.R {input}
		bash -e format_ucsc.sh
		"""
