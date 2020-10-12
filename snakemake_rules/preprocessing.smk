#### preprocessing rules ####

def get_value_from_df(wildcards, column_name):
	
	# get a value from the row of the df corresponding to this sample and dataset
	unique = f"{wildcards.dset}+++{wildcards.samp}"

	return toDo.loc[(toDo['unique'] == unique).idxmax(), column_name] 


rule check_bam_input_is_paired:
	input: 
		bam = lambda wildcards: get_value_from_df(wildcards, 'bam_file')
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
		bam = lambda wildcards: get_value_from_df(wildcards, 'bam_file'),
		ok = rules.check_bam_input_is_paired.output.ok
	output:
		r1 = temp("{outpath}/{dset}/reads/{samp}_1.fq.gz"),
		r2 = temp("{outpath}/{dset}/reads/{samp}_2.fq.gz"),
	conda:
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	shell:
		"""
		samtools view -b -F '0x900' {input.bam} |\
		samtools collate -O - |\
		samtools fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -
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
		r1 = lambda wildcards: get_for_seqprep(wildcards, '1'),
		r2 = lambda wildcards: get_for_seqprep(wildcards, '2')
	output:
		merged = "{outpath}/{dset}/merged_reads/{samp}.SeqPrep_merged.fastq.gz",
		proc_r1 = "{outpath}/{dset}/merged_reads/{samp}.1.fastq.gz",
		proc_r2 = "{outpath}/{dset}/merged_reads/{samp}.2.fastq.gz"
	conda:	
		"../envs/seqprep.yml"
	container:
		"docker://szsctt/seqprep:1"
	params:
		A = lambda wildcards: get_value_from_df(wildcards, "adapter_1"),
		B = lambda wildcards: get_value_from_df(wildcards, "adapter_2")
	shell:
		"""
		SeqPrep -A {params.A} -B {params.B} -f {input.r1} -r {input.r2} -1 {output.proc_r1} -2 {output.proc_r2} -s {output.merged}
		"""

rule touch_merged:
# if we don't want to do seqPrep, we still need to have an empty file of unmerged reads
	input:
		r1 = lambda wildcards: get_for_seqprep(wildcards, '1'),
		r2 = lambda wildcards: get_for_seqprep(wildcards, '2')
	output:
		merged = temp("{outpath}/{dset}/combined_reads/{samp}.mockMerged.fastq.gz")
	container:
		"docker://ubuntu:18.04"	
	shell:
		"""
		touch {output.merged}
		"""
		
		