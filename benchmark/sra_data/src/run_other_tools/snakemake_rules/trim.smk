rule trim:
	input:
		r1 = lambda wildcards: get_input_reads(wildcards, analysis_df, 1),
		r2 = lambda wildcards: get_input_reads(wildcards, analysis_df, 2)
	output:
		proc_r1 = temp("{outpath}/{dset}/trimmed_reads/{samp}.1.fastq.gz"),
		proc_r2 = temp("{outpath}/{dset}/trimmed_reads/{samp}.2.fastq.gz")
	resources:
		mem_mb= lambda wildcards, attempt: int(attempt * 5000),
		time = lambda wildcards, attempt: ('30:00', '2:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1
	conda:	
		"../envs/seqprep.yml"
	container:
		"docker://szsctt/seqprep:1"
	params:
		A = lambda wildcards: analysis_df_value(wildcards, analysis_df, 'adapter_1'),
		B = lambda wildcards: analysis_df_value(wildcards, analysis_df, 'adapter_2')
	shell:
		"""
		SeqPrep -A {params.A} -B {params.B} -f {input.r1} -r {input.r2} -1 {output.proc_r1} -2 {output.proc_r2}
		"""

