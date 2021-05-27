#### preprocessing rules ####
# rules should execute in the following order:
# 1. (if bam input) rule check_bam_input_is_paired: check if reads are paired, and number of r1 equals number of r2
# 2. (if bam input) rule bam_to_fastq: extract reads from bam file
# 3. (if dedup specified in config file) rule dedupe: remove duplicate pairs using clumpify dedupe
# 4. rule count_fastq: count the number of reads in r1 fastq to determine how to split files
# 5. rule split: split fastq files into smaller files to improve paralleization
# 6. rule SeqPrep or (rule SeqPrep_unmerged and rule touch_merged): depending on the trim and merge options specified in config: 
#		if trim only specified, then (rule SeqPrep_unmerged and rule touch_merged)
#		if trim and merge, or just merge specified then rule SeqPrep
#		if neither specified, then neither rule is performed

import functools
import re


# fuction to get mem_mb based on size of input files
def resources_list_with_min_and_max(file_name_list, attempt, mult_factor=2, minimum = 100, maximum = 120000):
	
	# get sum of size of files in file_name_list
	try:
		resource = int(sum([file.size for file in file_name_list]) / 1024 / 1024) * attempt * mult_factor
	# sometimes this doesn't work - not sure why...
	except WorkflowError:
		
		print(f"warning: couldn't get size of input files: using minimum {minimum}")
		resource = minimum * attempt

	resource = min(maximum, resource)
	
	return int(max(minimum, resource))
	
# function to get a value from the toDo df based on wildcards (to get correct row) and column_name
def get_value_from_df(wildcards, column_name):
	
	# get a value from the row of the df corresponding to this sample and dataset
	unique = f"{wildcards.dset}+++{wildcards.samp}"
	if column_name == 'align_cpus':
		return int(toDo.loc[(toDo['unique'] == unique).idxmax(), column_name])
		 
	return toDo.loc[(toDo['unique'] == unique).idxmax(), column_name] 
	
# number of cpus to use
cpus = functools.partial(get_value_from_df, column_name='align_cpus')
		
#get input fastq files depending on bam or fastq input
def get_reads_input(wildcards, read_type):

	assert read_type in ("1", "2")
	bam_suffix = get_value_from_df(wildcards, 'bam_file')

	# pass either reads extracted from bam or 
	if bam_suffix != "":
		if read_type == "1":
			return rules.bam_to_fastq.output.r1
		else:
			return rules.bam_to_fastq.output.r2
	else:
		if read_type == "1":
			return get_value_from_df(wildcards, 'R1_file')
		else:
			return get_value_from_df(wildcards, 'R2_file')

		
# Input functions for if had a bam or fastq as input
def get_for_split(wildcards, read_type):

	assert read_type in ("1", "2")
	bam_suffix = get_value_from_df(wildcards, 'bam_file')
	dedup = bool(get_value_from_df(wildcards, 'dedup'))

	# if we did de-duplication, then use the de-duplicated reads
	if dedup:
			if read_type == "1":
				return rules.dedupe.output.r1_dedup
			else:
				return rules.dedupe.output.r2_dedup
	
	# if no de-duplication, get input fastq files or reads extracted from bam
	else:
		return get_reads_input(wildcards, read_type)

# Get list of split values for sample
def get_split(wildcards):
	n_parts =  get_value_from_df(wildcards, "split")
	return range(int(n_parts))
	
# which cat to use in split rule
def get_cat(wildcards):
	dedup = bool(get_value_from_df(wildcards, 'dedup'))
	bam_suffix = get_value_from_df(wildcards, 'bam_file')
	
	# if we did de-duplication, reads will always be uncompressed
	if dedup:
		return 'cat'
	# same if we extracted reads from bam file
	elif bam_suffix != "":
		return 'cat'
	# otherwise, use original compression for input reads
	else:
		return get_value_from_df(wildcards, "cat")

def strip_wildcard_constraints(string_with_constraints):
	# when using the syntax rules.<rule>.output.<filename> to get the output of one rule for input to another
	# snakemake adds the wildcard constrains in after a comma
	# if we want to use these these as input to other rules, we need to strip out the commas
	
	string = string_with_constraints
	
	if re.search("\{\w+,\w+\}", string):
		string = re.sub(",.+?\}", "}", string)
	
	return string

rule write_analysis_summary:
	output:
		tsv = "{outpath}/summary/{dset}.analysis_conditions.tsv"
	run:
		toDo[toDo['dataset'] == wildcards.dset].to_csv(output.tsv, sep = "\t") 

rule check_bam_input_is_paired:
	input: 
		bam = lambda wildcards: get_value_from_df(wildcards, 'bam_file')
	output:
		ok = temp("{outpath}/{dset}/reads/{samp}.tmp"),
	resources:
		mem_mb=lambda wildcards, attempt, input: resources_list_with_min_and_max(input, attempt, 1.2, 300),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
		nodes = 1
	conda:
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	shell:
		"""
		FWD=$(samtools view -c -f 0x40 {input})
		if [[ "$FWD" == "0" ]]; then
			echo "Data must be paired-end"
			exit 1
		fi
		REV=$(samtools view -c -f 0x80 {input})
		if [[ "$FWD" != "$REV" ]]; then
			echo "Number of forward reads ($FWD) must match number of reverse reads ($REV)"
			exit 1
		fi
		touch {output.ok}
		"""	

rule bam_to_fastq:
	input:
		bam = lambda wildcards: get_value_from_df(wildcards, 'bam_file'),
		ok = rules.check_bam_input_is_paired.output.ok
	output:
		r1 = temp("{outpath}/{dset}/reads/{samp}_1.fq"),
		r2 = temp("{outpath}/{dset}/reads/{samp}_2.fq"),
	resources:
		mem_mb=lambda wildcards, attempt, input: resources_list_with_min_and_max(input, attempt, 1.2, 1000),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
		nodes = 1
	conda:
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	resources:
		mem_mb=lambda wildcards, attempt, input: resources_list_with_min_and_max(input, attempt)
	shell:
		"""
		samtools view -b -F '0x900' {input.bam} |\
		samtools collate -O - |\
		samtools fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -
		"""
		
rule dedupe:
	input:
		r1 = lambda wildcards: get_reads_input(wildcards, "1"),
		r2 = lambda wildcards: get_reads_input(wildcards, "2")
	output:
		r1_dedup = temp("{outpath}/{dset}/dedup_reads/{samp}_1.fq"),
		r2_dedup = temp("{outpath}/{dset}/dedup_reads/{samp}_2.fq")
	params:
		n_subs = lambda wildcards: get_value_from_df(wildcards, "dedup_subs"),
		mem_mb = lambda wildcards, resources: max(int(resources.mem_mb * 0.8), 1000)
	threads: cpus
	conda:	
		"../envs/bbmap.yml"
	container:
		"docker://szsctt/bbmap:1"	
	resources:
		mem_mb = lambda wildcards, attempt, input: resources_list_with_min_and_max(input, attempt, 8),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
		nodes = 1 # need this for pearcey so that job doesn't get split over multiple nodes
	shell:
		"""
		clumpify.sh -Xmx{params.mem_mb}m in1={input.r1} in2={input.r2} out1={output.r1_dedup} out2={output.r2_dedup} dedupe=t ac=f subs={params.n_subs}{threads}
		"""

rule count_fastq:
	input:
		reads = lambda wildcards: get_for_split(wildcards, '1')
	output:
		count_reads = temp("{outpath}/{dset}/reads/{samp}_count.tmp")	
	params:
		n_total =  lambda wildcards: get_value_from_df(wildcards, "split"),
		cat = lambda wildcards: get_cat(wildcards),
	resources:
		mem_mb = lambda wildcards, attempt, input: resources_list_with_min_and_max(input, attempt, 0.5, 500, 5000),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
		nodes = 1
	shell:
		"""
		rm -f {output.count_reads}
		count=$({params.cat} {input.reads} | wc -l)
		split_n={params.n_total}
		chunk_lines=$((count / split_n))
		curr_start=1
		while [ $curr_start -le $count ]
		do
			end=$(( 4 * ((curr_start-1)+chunk_lines) / 4))
			mod=$((end % 4))
			end=$((end + 4 - mod))
			echo $curr_start','$end >> {output.count_reads}
			curr_start=$((end + 1))
		done
		"""

rule split_fastq:
	input:
		reads = lambda wildcards: get_for_split(wildcards, wildcards.read_num),
		count_reads = "{outpath}/{dset}/reads/{samp}_count.tmp"	
	output:
		reads = temp("{outpath}/{dset}/split_reads/{samp}_{read_num}.{part}.fq"),
	params:
		n_total =  lambda wildcards: get_value_from_df(wildcards, "split"),
		cat = lambda wildcards: get_cat(wildcards),
	resources:
		mem_mb = lambda wildcards, attempt, input: resources_list_with_min_and_max(input, attempt, 0.5, 500, 5000),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
		nodes = 1
	wildcard_constraints:
		read_num = "1|2",
		part = "\d+"
	shell:
		"""
		if [[ {params.n_total} -eq 1 ]]
		then
			{params.cat} {input.reads} > {output.reads}
		else
			part=$(({wildcards.part}+1))
			lines=$(sed -n $part" p" {input.count_reads})
			cmd="{params.cat} {input.reads} | sed -n  '"$lines" p' > {output.reads}"
			eval $cmd
		fi
		"""

rule seqPrep:
# if we're doing it
	input:
		r1 = "{outpath}/{dset}/split_reads/{samp}_1.{part}.fq",
		r2 = "{outpath}/{dset}/split_reads/{samp}_2.{part}.fq"
	output:
		merged = temp("{outpath}/{dset}/merged_reads/{samp}.{part}.SeqPrep_merged.fastq.gz"),
		proc_r1 = temp("{outpath}/{dset}/merged_reads/{samp}.{part}.1.fastq.gz"),
		proc_r2 = temp("{outpath}/{dset}/merged_reads/{samp}.{part}.2.fastq.gz")
	conda:	
		"../envs/seqprep.yml"
	container:
		"docker://szsctt/seqprep:1"
	resources:
		mem_mb = lambda wildcards, attempt, input: resources_list_with_min_and_max(input, attempt, 0.5, 500, 5000),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
		nodes = 1
	params:
		A = lambda wildcards: get_value_from_df(wildcards, "adapter_1"),
		B = lambda wildcards: get_value_from_df(wildcards, "adapter_2")
	shell:
		"""
		SeqPrep -A {params.A} -B {params.B} -f {input.r1} -r {input.r2} -1 {output.proc_r1} -2 {output.proc_r2} -s {output.merged}
		"""
		
rule seqPrep_unmerged:
	input:
		r1 = "{outpath}/{dset}/split_reads/{samp}_1.{part}.fq",
		r2 = "{outpath}/{dset}/split_reads/{samp}_2.{part}.fq"
	output:
		proc_r1 = temp("{outpath}/{dset}/trimmed_reads/{samp}.{part}.1.fastq.gz"),
		proc_r2 = temp("{outpath}/{dset}/trimmed_reads/{samp}.{part}.2.fastq.gz")
	conda:	
		"../envs/seqprep.yml"
	container:
		"docker://szsctt/seqprep:1"
	resources:
		mem_mb = lambda wildcards, attempt, input: resources_list_with_min_and_max(input, attempt, 0.5, 500, 5000),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
		nodes = 1
	params:
		A = lambda wildcards: get_value_from_df(wildcards, "adapter_1"),
		B = lambda wildcards: get_value_from_df(wildcards, "adapter_2")
	shell:
		"""
		SeqPrep -A {params.A} -B {params.B} -f {input.r1} -r {input.r2} -1 {output.proc_r1} -2 {output.proc_r2}
		"""

# if we don't want to do merging, we still need to have an empty file of unmerged reads
rule touch_merged:
	output:
		merged = temp("{outpath}/{dset}/combined_reads/{samp}.{part}.mockMerged.fastq.gz")
	container:
		"docker://ubuntu:18.04"	
	shell:
		"""
		touch {output.merged}
		"""
		
		
