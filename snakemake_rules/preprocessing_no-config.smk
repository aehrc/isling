#### preprocessing rules ####
import pdb

def get_data(wildcards, column_name):

	row = toDo.loc[toDo['unique'] == f"{wildcards.dset}+++{wildcards.samp}"]
		
	return row[column_name][0]

rule check_bam_input_is_paired:
	input: lambda wildcards: get_data(wildcards, "bam_file")
	output:
		ok = temp("{outpath}/{dset}/reads/{samp}.tmp"),
	conda:
		"envs/bwa.yml"
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
		bam = lambda wildcards: get_data(wildcards, "bam_file"),
		ok = lambda wildcards: f"{wildcards.outpath}/{wildcards.dset}/reads/{wildcards.samp}.tmp"
	output:
		sorted_bam = temp("{outpath}/{dset}/reads/{samp}.sorted.bam"),
		r1 = temp("{outpath}/{dset}/reads/{samp}_1.fq.gz"),
		r2 = temp("{outpath}/{dset}/reads/{samp}_2.fq.gz"),
	conda:
		"envs/picard.yml"
	container:
		"docker://szsctt/picard:1"
	shell:
		"""
		picard SortSam I={input.bam} O={output.sorted_bam} SORT_ORDER=queryname
		picard SamToFastq I={output.sorted_bam} FASTQ={output.r1} SECOND_END_FASTQ={output.r2}
		"""
		
def get_for_merge(wildcards, column):

	row = toDo.loc[toDo['unique'] == f"{wildcards.dset}+++{wildcards.samp}"]
	bam = bool(row['bam_file'][0]) != ""
	
	# if we started from a bam file, return the output of the conversion rule
	if bam:
		if column == "R1_file":
			return rules.bam_to_fastq.output.r1
		else:
			return rules.bam_to_fastq.output.r2
	# if we started from reads
	else:
		return row[column][0]
	

rule seqPrep:
# if we're doing it
	input:
		r1 = lambda wildcards: get_for_merge(wildcards, 'R1_file'),
		r2 = lambda wildcards: get_for_merge(wildcards, 'R2_file')
	output:
		merged = temp("{outpath}/{dset}/merged_reads/{samp}.SeqPrep_merged.fastq.gz"),
		proc_r1 = temp("{outpath}/{dset}/merged_reads/{samp}.1.fastq.gz"),
		proc_r2 = temp("{outpath}/{dset}/merged_reads/{samp}.2.fastq.gz"),
		all = temp("{outpath}/{dset}/merged_reads/{samp}.all.fastq.gz")
	conda:	
		"envs/seqprep.yml"
	container:
		"docker://szsctt/seqprep:1"
	params:
		A = lambda wildcards: f"-A {get_data(wildcards, 'adapter_1')}",
		B = lambda wildcards: f"-B {get_data(wildcards, 'adapter_2')}"
	shell:
		"""
		SeqPrep {params} -f {input.r1} -r {input.r2} -1 {output.proc_r1} -2 {output.proc_r2} -s {output.merged}
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
		r1 = lambda wildcards: get_for_merge(wildcards, 'R1_file'),
		r2 = lambda wildcards: get_for_merge(wildcards, 'R2_file')
	output:
		proc_r1 = temp("{outpath}/{dset}/combined_reads/{samp}.1.fastq.gz"),
		proc_r2 = temp("{outpath}/{dset}/combined_reads/{samp}.2.fastq.gz"),
		all = temp("{outpath}/{dset}/combined_reads/{samp}.all.fastq.gz")
	params:
		cat = lambda wildcards, input: get_compression(input),
	container:
		"docker://ubuntu:18.04"	
	shell:
		"""
		{params.cat} {input.r1} | gzip > {output.proc_r1}
		{params.cat} {input.r2} | gzip > {output.proc_r2}
		{params.cat} {input.r1} {input.r2} | gzip > {output.all}
		"""
