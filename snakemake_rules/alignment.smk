

#functions for if we did seqPrep or not
def get_for_align(wildcards, read_type):

	assert read_type in {'unmerged_r1', 'unmerged_r2', 'merged'}

	merge = bool(get_value_from_df(wildcards, 'merge'))
	trim = bool(get_value_from_df(wildcards, 'trim'))
	dedup = bool(get_value_from_df(wildcards, 'dedup'))
	
	# did we merge R1 and R2?
	if merge is True:
		if read_type == 'unmerged_r1':
			return rules.seqPrep.output.proc_r1
		if read_type == 'unmerged_r2':
			return rules.seqPrep.output.proc_r2
		else:
			return rules.seqPrep.output.merged
	# if we trimmed but didn't merge
	elif trim is True:
		if read_type == 'unmerged_r1':
			return rules.seqPrep_unmerged.output.proc_r1
		if read_type == 'unmerged_r2':
			return rules.seqPrep_unmerged.output.proc_r2
		else:
			return rules.touch_merged.output.merged
	# if we did de-duplication, but didn't trim or merge
	elif dedup is True:
		if read_type == 'unmerged_r1':
			return rules.dedupe.output.r1_dedup
		if read_type == 'unmerged_r2':
			return rules.dedupe.output.r2_dedup
		else:
			return rules.touch_merged.output.merged
	# if we didn't do any of these
	else:
		if read_type == 'unmerged_r1':
			return "{outpath}/{dset}/split_reads/{samp}_1.{part}.fq"
		if read_type == 'unmerged_r2':
			return "{outpath}/{dset}/split_reads/{samp}_2.{part}.fq"
		else:
			return rules.touch_merged.output.merged


#### alignments ####
rule index:
	input:
		fa = lambda wildcards: ref_names[wildcards.genome]
	output:
		temp(expand("{outpath}/references/{genome}/{genome}.{ext}", ext=["ann", "amb", "bwt", "pac", "sa"], allow_missing=True))
	conda: 
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	params:
		prefix = lambda wildcards, output: path.splitext(output[0])[0]
	resources:
		mem_mb=lambda wildcards, attempt, input: resources_list_with_min_and_max((input.fa,), attempt, 5, 1000)
	shell:
		"bwa index -p {params.prefix} {input.fa}"


rule align_bwa_virus:
	input:
		idx = lambda wildcards: multiext(get_value_from_df(wildcards, 'virus_prefix'), ".ann", ".amb", ".bwt", ".pac", ".sa"),
		merged = lambda wildcards: get_for_align(wildcards, "merged"),
		r1 = lambda wildcards: get_for_align(wildcards, "unmerged_r1"),
		r2 = lambda wildcards: get_for_align(wildcards, "unmerged_r2"),
	output:
		single = temp("{outpath}/{dset}/virus_aligned/{samp}.{part}.{virus}.bwaSingle.sam"),
		paired = temp("{outpath}/{dset}/virus_aligned/{samp}.{part}.{virus}.bwaPaired.sam"),
		combined = temp("{outpath}/{dset}/virus_aligned/{samp}.{part}.{virus}.sam"),
	params:
		index = lambda wildcards, input: os.path.splitext(input.idx[0])[0],
		mapping = lambda wildcards: get_value_from_df(wildcards, 'bwa_mem_params'),
		single_RG = lambda wildcards: f"-R '@RG\\tID:{wildcards.samp}_{wildcards.virus}_merged\\tSM:{wildcards.samp}\\tPM:merged'",
		paired_RG = lambda wildcards: f"-R '@RG\\tID:{wildcards.samp}_{wildcards.virus}_unmerged\\tSM:{wildcards.samp}\\tPM:unmerged'"
	resources:
		mem_mb=lambda wildcards, attempt, input: resources_list_with_min_and_max(input.idx, attempt, 5, 2000)
	conda:
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	threads: cpus
	shell:
		"""
		bwa mem -t {threads} {params.mapping} {params.single_RG} -o {output.single} {params.index} {input.merged}
		
		bwa mem -t {threads} {params.mapping} {params.paired_RG} -o {output.paired} {params.index} {input.r1} {input.r2} 
		
		samtools merge {output.combined} {output.single} {output.paired}
		"""

def get_sam(wildcards, readType, genome):
	
	assert readType in ['single', 'paired', 'combined']
	assert genome in ['host', 'virus']

	merge = bool(get_value_from_df(wildcards, 'merge'))
	
	# if we want host alignment
	if genome == "virus":
		if readType == "single":
			return rules.align_bwa_virus.output.single
		elif readType == "paired":
			return rules.align_bwa_virus.output.paired
		else:
			return rules.align_bwa_virus.output.combined
	# if we want the host alignment
	else:
		if readType == "single":
			return rules.align_bwa_host_single.output.sam
		elif readType == "paired":
			return rules.align_bwa_host_paired.output.sam
		else:
			return rules.combine_host.output.combined
		
rule extract_to_fastq_single:
	input:
		aligned = lambda wildcards: get_sam(wildcards, "single", "virus"),
	output:
		fastq = temp("{outpath}/{dset}/virus_aligned/{samp}.{part}.bwaSingle.mappedTo{virus}.fastq.gz"),
	group: "host_align"
	conda:
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	shell:
		"""
		# 0x4 - read unmapped
		# 0x100 - not primary alignment
		# 0x800 - secondary alignment
		samtools view -h -F 0x4 -F 0x800 -F 0x100 -o - {input.aligned} | samtools fastq -0 {output.fastq}  
		"""

rule extract_vAligned_paired:
	input:
		aligned = lambda wildcards: get_sam(wildcards, "paired", "virus")
	output:
		pvBam_readMap_mateUnmap = temp("{outpath}/{dset}/virus_aligned/{samp}.{part}.{virus}.bwaPaired.L.bam"),
		pvBam_readUnmap_mateMap = temp("{outpath}/{dset}/virus_aligned/{samp}.{part}.{virus}.bwaPaired.R.bam"),
		pvBam_bothMapped = temp("{outpath}/{dset}/virus_aligned/{samp}.{part}.{virus}.bwaPaired.B.bam"),
		bam = temp("{outpath}/{dset}/virus_aligned/{samp}.{part}.{virus}.bwaPaired.mapped.bam"),
	group: "host_align"
	conda:
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	shell:
		"""
		samtools view -hb -F 0x4 -f 0x8 -F 0x800 -o {output.pvBam_readMap_mateUnmap} {input.aligned}
		samtools view -hb -f 0x4 -F 0x8 -F 0x800 -o {output.pvBam_readUnmap_mateMap} {input.aligned}
		samtools view -hb -F 0x4 -F 0x8 -F 0x800 -o {output.pvBam_bothMapped} {input.aligned}
		samtools merge {output.bam} {output.pvBam_readMap_mateUnmap} {output.pvBam_bothMapped} {output.pvBam_readUnmap_mateMap}
		"""
		
rule extract_to_fastq_paired:
	input:
		bam = rules.extract_vAligned_paired.output.bam
	output:
		fastq1 = temp("{outpath}/{dset}/virus_aligned/{samp}.{part}.bwaPaired.mappedTo{virus}.1.fastq.gz"),
		fastq2 = temp("{outpath}/{dset}/virus_aligned/{samp}.{part}.bwaPaired.mappedTo{virus}.2.fastq.gz")
	group: "host_align"
	conda:
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	shell:
		"""
		samtools collate -O {input.bam} | samtools fastq -1 {output.fastq1} -2 {output.fastq2}
		"""

rule align_bwa_host_single:
	input:	
		idx = lambda wildcards: multiext(get_value_from_df(wildcards, 'host_prefix'), ".ann", ".amb", ".bwt", ".pac", ".sa"),
		all = rules.extract_to_fastq_single.output.fastq,
	output:
		sam = temp("{outpath}/{dset}/host_aligned/{samp}.{part}.{host}.readsFrom{virus}.bwaSingle.sam"),
	group: "host_align"
	conda: 
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	resources:
		mem_mb=lambda wildcards, attempt, input: resources_list_with_min_and_max(input.idx, attempt, 5, 2000)
	params:
		index = lambda wildcards, input: os.path.splitext(input.idx[0])[0],
		mapping = lambda wildcards: get_value_from_df(wildcards, 'bwa_mem_params'),
		RG = lambda wildcards: f"-R '@RG\\tID:{wildcards.samp}_{wildcards.host}_merged\\tSM:{wildcards.samp}\\tPM:merged'"
	threads: cpus
	shell:		
		"""
		bwa mem -t {threads} {params.mapping} {params.RG} -o {output.sam} {params.index} {input.all}
		"""
		
rule align_bwa_host_paired:
	input:	
		idx = lambda wildcards: multiext(get_value_from_df(wildcards, 'host_prefix'), ".ann", ".amb", ".bwt", ".pac", ".sa"),
		r1 = rules.extract_to_fastq_paired.output[0],
		r2 = rules.extract_to_fastq_paired.output[1]
	output:
		sam = temp("{outpath}/{dset}/host_aligned/{samp}.{part}.{host}.readsFrom{virus}.bwaPaired.sam"),
	group: "host_align"
	conda: 
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	resources:
		mem_mb=lambda wildcards, attempt, input: resources_list_with_min_and_max(input.idx, attempt, 5, 2000)
	params:
		index = lambda wildcards, input: os.path.splitext(input.idx[0])[0],
		mapping = lambda wildcards: get_value_from_df(wildcards, 'bwa_mem_params'),
		RG = lambda wildcards: f"-R '@RG\\tID:{wildcards.samp}_{wildcards.host}_unmerged\\tSM:{wildcards.samp}\\tPM:unmerged'"
	threads: cpus
	shell:		
		"""
		bwa mem -t {threads} {params.mapping} {params.RG} -o {output.sam} {params.index} {input.r1} {input.r2}
		"""
		
rule combine_host:
	input:
		paired = rules.align_bwa_host_paired.output.sam,
		single = rules.align_bwa_host_single.output.sam
	output:
		combined = temp("{outpath}/{dset}/host_aligned/{samp}.{part}.{host}.readsFrom{virus}.sam")
	group: "host_align"
	conda: 
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	resources:
		mem_mb=lambda wildcards, attempt, input: resources_list_with_min_and_max((input.paired, input.single), attempt)
	shell:		
		"""
		samtools merge {output.combined} {input.single} {input.paired}
		"""
	

#### sam file manipulations ####
rule merge_virus_sams:
	message: "Merging virus sam files to one single sam file."
	input:
		lambda wildcards: expand("{{outpath}}/{{dset}}/virus_aligned/{{samp}}.{parts}.{{virus}}.sam", parts = get_split(wildcards))
	output:
		temp("{outpath}/{dset}/virus_aligned/{samp}.{virus}.sam")
	conda:
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	threads: 1
	shell:
		"""
		samtools merge {output} {input}
		"""

rule merge_host_sams:
	message: "Merging host sam files to one single sam file."
	input:
		lambda wildcards: expand("{{outpath}}/{{dset}}/host_aligned/{{samp}}.{parts}.{{host}}.readsFrom{{virus}}.bwaSingle.sam", parts = get_split(wildcards))
	output:
		temp("{outpath}/{dset}/host_aligned/{samp}.{host}.readsFrom{virus}.bwaSingle.sam")
	conda:
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	threads: 1
	shell:
		"""
		samtools merge {output} {input}
		"""

rule convert_virus_sam_to_bam:
	input:
		sam = rules.merge_virus_sams.output
	output:
		bam = "{outpath}/{dset}/virus_aligned/{samp}.{virus}.bam",
		bai = "{outpath}/{dset}/virus_aligned/{samp}.{virus}.bam.bai"
	conda: 
		"../envs/bwa.yml"	
	container:
		"docker://szsctt/bwa:1"
#	resources: # removed due to error in splitFASTQ-feature
#		mem_mb=lambda wildcards, attempt, input: resources_list_with_min_and_max((input.sam,), attempt, 2, 1000, 10000)
	shell:
		"""
		samtools sort -o {output.bam} {input.sam}
		samtools index {output.bam}
		"""

rule convert_host_sam_to_bam:
	input:
		sam = rules.merge_host_sams.output
	output:
		bam = "{outpath}/{dset}/host_aligned/{samp}.{host}.readsFrom{virus}.bam",
		bai = "{outpath}/{dset}/host_aligned/{samp}.{host}.readsFrom{virus}.bam.bai"
	conda: 
		"../envs/bwa.yml"	
	container:
		"docker://szsctt/bwa:1"
	shell:
		"""
		samtools sort -o {output.bam} {input.sam}
		samtools index {output.bam}
		"""

