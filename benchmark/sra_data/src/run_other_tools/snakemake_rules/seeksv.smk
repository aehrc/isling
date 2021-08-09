
rule host_virus_index_seeksv:
	input:
		virus = lambda wildcards: ref_names[wildcards.virus],
		host = lambda wildcards: ref_names[wildcards.host]
	output:
		fa = "{outpath}/seeksv_refs/{virus}/{host}_{virus}.fas",
		idx = multiext("{outpath}/seeksv_refs/{virus}/{host}_{virus}.fas", ".amb", ".ann", ".bwt", ".pac", ".sa")
	container:
		"docker://szsctt/seeksv:1"
	resources:
		mem_mb= lambda wildcards, attempt, input: resources_list_with_min_and_max((input.host, input.virus), attempt, 5, 1000),
		time = lambda wildcards, attempt: ('2:00:00', '24:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1
	shell:
		"""
		cat {input.virus} {input.host} > {output.fa}
		bwa index {output.fa}
		"""
		
rule align_seeksv_all:
	input:
		fastq1 = lambda wildcards: get_reads(wildcards, analysis_df, rules, 1),
		fastq2 = lambda wildcards: get_reads(wildcards, analysis_df, rules, 2),
		idx = rules.host_virus_index_seeksv.output.idx
	output:
		bam = "{outpath}/{dset}/seeksv/aln/{samp}.{host}.{virus}.bam",
		bai = "{outpath}/{dset}/seeksv/aln/{samp}.{host}.{virus}.bam.bai"
	params:
		prefix = lambda wildcards, input: os.path.splitext(input.idx[0])[0]	
	threads: 8
	container:
		"docker://szsctt/seeksv:1"
	resources:
		mem_mb= lambda wildcards, attempt, input: resources_list_with_min_and_max(input.idx, attempt, 5, 1000),
		nodes = 1,
		time = lambda wildcards, attempt: ('2:00:00', '24:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
	shell:
		"""
		rm -f {output.bam}.tmp*
		bwa mem -t {threads} {params.prefix} {input.fastq1} {input.fastq2} | samtools sort -o {output.bam} -
		samtools index {output.bam}
		"""
		
rule dedup_seeksv:
	input:
		bam = rules.align_seeksv_all.output.bam
	output:
		bam = "{outpath}/{dset}/seeksv/aln/{samp}.{host}.{virus}.dup.bam",
		metrics = "{outpath}/{dset}/seeksv/aln/{samp}.{host}.{virus}.dup.txt"
	resources:
		mem_mb= lambda wildcards, attempt, input: resources_list_with_min_and_max((input.bam, ), attempt, mult_factor=1, minimum=3000),
		time = lambda wildcards, attempt: ('30:00', '2:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1
	params:
		mem_gb_sort = lambda wildcards, resources: int(resources.mem_mb / 1e3 / 3),
		mem_gb_dup = lambda wildcards, resources: int(resources.mem_mb / 1e3 / 3),
	container:
		"docker://szsctt/seeksv:1"	
	shell:
		"""
		java -Xmx{params.mem_gb_sort}g -jar ${{PICARD}} SortSam \
			I={input.bam} \
			VALIDATION_STRINGENCY=LENIENT \
			COMPRESSION_LEVEL=0 \
			O=/dev/stdout \
			SORT_ORDER=queryname | java -Xmx{params.mem_gb_dup}g -jar ${{PICARD}} MarkDuplicates \
			I=/dev/stdin \
			VALIDATION_STRINGENCY=LENIENT \
			METRICS_FILE={output.metrics} \
			O={output.bam}
		"""
		
def get_seeksv_alignment(wildcards):
	# if we want to do dedupliation
	if analysis_df_tool_value(wildcards, analysis_df, 'seeksv', 'dedup') == 1:
		return rules.dedup_seeksv.output.bam
	
	# no deduplication
	return rules.align_seeksv_all.output.bam
		
rule sort_seeksv:
	input:
		bam = lambda wildcards: get_seeksv_alignment(wildcards)
	output:
		bam = "{outpath}/{dset}/seeksv/aln/{samp}.{host}.{virus}.sorted.bam",
		bai = "{outpath}/{dset}/seeksv/aln/{samp}.{host}.{virus}.sorted.bam.bai"
	resources:
		mem_mb= lambda wildcards, attempt, input: resources_list_with_min_and_max((input.bam, ), attempt, minimum=1000),
		time = lambda wildcards, attempt: ('30:00', '2:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1
	params:
		mem_gb = lambda wildcards, resources: int(resources.mem_mb / 1e3)
	container:
		"docker://szsctt/seeksv:1"	
	shell:
		"""
		java -Xmx{params.mem_gb}g -jar ${{PICARD}} SortSam \
			I={input.bam} \
			VALIDATION_STRINGENCY=LENIENT \
			COMPRESSION_LEVEL=0 \
			O={output.bam} \
			SORT_ORDER=coordinate
			
		samtools index {output.bam}
		"""
		
rule seeksv_getclip:
	input:
		bam = rules.sort_seeksv.output.bam,
		bai = rules.sort_seeksv.output.bai
	output: 
		clip_fq = "{outpath}/{dset}/seeksv/clipped_reads/{samp}.{host}.{virus}.clip.fq.gz",
		clip = "{outpath}/{dset}/seeksv/clipped_reads/{samp}.{host}.{virus}.clip.gz",
		unmapped_1 = "{outpath}/{dset}/seeksv/clipped_reads/{samp}.{host}.{virus}.unmapped_1.fq.gz",
		unmapped_2 = "{outpath}/{dset}/seeksv/clipped_reads/{samp}.{host}.{virus}.unmapped_2.fq.gz",	
	resources:
		mem_mb= lambda wildcards, attempt, input: resources_list_with_min_and_max((input.bam, ), attempt, 1),
		time = lambda wildcards, attempt: ('30:00', '2:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1
	params:
		prefix = lambda wildcards, output: os.path.splitext(os.path.splitext(output.clip)[0])[0]
	container:
		"docker://szsctt/seeksv:1"
	shell:
		"""
		/var/work/seeksv/seeksv getclip -o {params.prefix} {input.bam}
		"""
		
rule align_seeksv_clip:
	input:
		fq = rules.seeksv_getclip.output.clip_fq,
		idx = rules.host_virus_index_seeksv.output.idx
	output:
		bam = "{outpath}/{dset}/seeksv/aln/{samp}.{host}.{virus}.clip.bam"
	resources:
		mem_mb= lambda wildcards, attempt, input: resources_list_with_min_and_max(input.idx, attempt, 5, 1000),
		time = lambda wildcards, attempt: ('30:00', '2:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1
	threads: 8
	params:
		prefix = lambda wildcards, input: os.path.splitext(input.idx[0])[0]
	container:
		"docker://szsctt/seeksv:1"
	shell:
		"""
		bwa mem -t {threads}  {params.prefix} {input.fq} | samtools view  -Sb -o {output.bam} -
		"""

rule seeksv:
	input:
		clip = rules.seeksv_getclip.output.clip,
		bam = rules.align_seeksv_all.output.bam,
		bam_clip = rules.align_seeksv_clip.output.bam
	output:
		ints = "{outpath}/{dset}/seeksv/ints/{samp}.{host}.{virus}.integrations.txt",
		unmapped = "{outpath}/{dset}/seeksv/ints/{samp}.{host}.{virus}.unmapped.clip.fq.gz"
	resources:
		mem_mb= lambda wildcards, attempt, input: resources_list_with_min_and_max((input.clip, input.bam, input.bam_clip), attempt),
		time = lambda wildcards, attempt: ('30:00', '2:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1
	container:
		"docker://szsctt/seeksv:1"
	shell:
		"""
		/var/work/seeksv/seeksv getsv {input.bam_clip} {input.bam} {input.clip} {output.ints} {output.unmapped}
		"""		

