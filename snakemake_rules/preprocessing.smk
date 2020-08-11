#### preprocessing rules ####

def get_value_from_df(wildcards, column_name):
	
	# get a value from the row of the df corresponding to this sample and dataset
	unqiue = f"{wildcards.dset}+++{wildcards.samp}"
	return toDo.loc[(toDo['unique'] == unique).idxmax(), column_name]


rule check_bam_input_is_paired:
	input: lambda wildcards: get_value_from_df(wildcards, 'bam_file')
	output:
		ok = temp("{outpath}/{dset}/reads/{samp}.tmp"),
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
			echo "Number of forward reads must match number of paired reads"
			exit 1
		fi
		touch {output.ok}
		"""	

rule bam_to_fastq:
	input:
		bam = lambda wildcards:  get_value_from_df(wildcards, 'bam_file'),
		ok = lambda wildcards: f"{wildcards.outpath}/{wildcards.dset}/reads/{wildcards.samp}.tmp"
	output:
		sorted_bam = temp("{outpath}/{dset}/reads/{samp}.sorted.bam"),
		r1 = temp("{outpath}/{dset}/reads/{samp}_1.fq.gz"),
		r2 = temp("{outpath}/{dset}/reads/{samp}_2.fq.gz"),
	conda:
		"../envs/picard.yml"
	container:
		"docker://szsctt/picard:1"
	shell:
		"""
		picard SortSam I={input.bam} O={output.sorted_bam} SORT_ORDER=queryname
		picard SamToFastq I={output.sorted_bam} FASTQ={output.r1} SECOND_END_FASTQ={output.r2}
		"""


# input functions for if had a bam or fastq as input
def get_for_seqprep(wildcards, read_type):

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


rule seqPrep:
# if we're doing it
	input:
		r1 = lambda wildcards: get_value_from_df(wildcards, 'R1_file'),
		r2 = lambda wildcards: get_value_from_df(wildcards, 'R2_file')
	output:
		merged = "{outpath}/{dset}/merged_reads/{samp}.SeqPrep_merged.fastq.gz",
		proc_r1 = "{outpath}/{dset}/merged_reads/{samp}.1.fastq.gz",
		proc_r2 = "{outpath}/{dset}/merged_reads/{samp}.2.fastq.gz"
	conda:	
		"../envs/seqprep.yml"
	container:
		"docker://szsctt/seqprep:1"
	params:
		A = lambda wildcards: config[wildcards.dset]["read1-adapt"],
		B = lambda wildcards: config[wildcards.dset]["read2-adapt"]
	shell:
		"""
		SeqPrep -A {params.A} -B {params.B} -f {input.r1} -r {input.r2} -1 {output.proc_r1} -2 {output.proc_r2} -s {output.merged}
		"""

rule touch_merged:
# if we don't want to do seqPrep, we still need to have an empty file of unmerged reads
	input:
		r1 = lambda wildcards: get_value_from_df(wildcards, 'R1_file'),
		r2 = lambda wildcards: get_value_from_df(wildcards, 'R2_file')
	output:
		merged = temp("{outpath}/{dset}/combined_reads/{samp}.mockMerged.fastq.gz")
	container:
		"docker://ubuntu:18.04"	
	shell:
		"""
		touch {output.merged}
		"""
		
		
