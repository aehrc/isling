#### preprocessing rules ####


rule check_bam_input_is_paired:
	input: lambda wildcards: path.normpath(f"{config[wildcards.dset]['read_folder']}/{wildcards.samp}{config[wildcards.dset]['bam_suffix']}")
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
		bam = lambda wildcards: path.normpath(f"{config[wildcards.dset]['read_folder']}/{wildcards.samp}{config[wildcards.dset]['bam_suffix']}"),
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
	row_idx = list(toDo.loc[:,'unique']).index(f"{wildcards.dset}+++{wildcards.samp}")
	# get true/True values for dedup
	if 'bam_suffix' in config[wildcards.dset]:
		if read_type == "1":
			return "{outpath}/{dset}/reads/{samp}_1.fq.gz"
		else:
			return "{outpath}/{dset}/reads/{samp}_2.fq.gz"
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
		r1 = lambda wildcards: get_for_seqprep(wildcards, "1"),
		r2 = lambda wildcards: get_for_seqprep(wildcards, "2")
	output:
		merged = temp("{outpath}/{dset}/combined_reads/{samp}.mockMerged.fastq.gz")
	container:
		"docker://ubuntu:18.04"	
	shell:
		"""
		touch {output.merged}
		"""
		
		
