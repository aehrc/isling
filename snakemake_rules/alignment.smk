
#functions for if we did seqPrep or not
def get_for_align(wildcards, read_type):

	assert read_type in {'unmerged_r1', 'unmerged_r2', 'merged'}

	merge = bool(get_value_from_df(wildcards, 'merge'))
	trim = bool(get_value_from_df(wildcards, 'trim'))
	
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
	# if we didn't do either
	else:
		if read_type == 'unmerged_r1':
			return get_for_seqprep(wildcards, "1")
		if read_type == 'unmerged_r2':
			return get_for_seqprep(wildcards, "2")
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
		mem_mb=lambda wildcards, attempt, input: attempt * 3 * int(os.stat(input.fa).st_size/1e6)
	shell:
		"bwa index -p {params.prefix} {input.fa}"


rule align_bwa_virus:
	input:
		idx = expand("{outpath}/references/{virus}/{virus}.{ext}", ext=["ann", "amb", "bwt", "pac", "sa"], allow_missing=True),
		merged = lambda wildcards: get_for_align(wildcards, "merged"),
		r1 = lambda wildcards: get_for_align(wildcards, "unmerged_r1"),
		r2 = lambda wildcards: get_for_align(wildcards, "unmerged_r2"),
	
	output:
		single = temp("{outpath}/{dset}/virus_aligned/{samp}.{virus}.bwaSingle.sam"),
		paired = temp("{outpath}/{dset}/virus_aligned/{samp}.{virus}.bwaPaired.sam"),
		combined = temp("{outpath}/{dset}/virus_aligned/{samp}.{virus}.sam"),
	params:
		index = lambda wildcards, input: path.splitext(input.idx[0])[0],
		mapping = lambda wildcards: get_value_from_df(wildcards, 'bwa_mem_params'),
		single_RG = lambda wildcards: f"-R '@RG\\tID:{wildcards.samp}_{wildcards.virus}_merged\\tSM:{wildcards.samp}\\tPM:merged'",
		paired_RG = lambda wildcards: f"-R '@RG\\tID:{wildcards.samp}_{wildcards.virus}_unmerged\\tSM:{wildcards.samp}\\tPM:unmerged'"
	resources:
		mem_mb=lambda wildcards, attempt, input: attempt * 3 * sum([int(os.stat(file).st_size/1e6) for file in input.idx])
	conda:
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	threads: 5
	shell:
		"""
		bwa mem -t {threads} {params.mapping} {params.single_RG} -o {output.single} {params.index} {input.merged}
		
		bwa mem -t {threads} {params.mapping} {params.paired_RG} -o {output.paired} {params.index} {input.r1} {input.r2} 
		
		samtools merge {output.combined} {output.single} {output.paired}
		"""
		

def get_sam(wildcards, readType, genome):

	#pdb.set_trace()
	
	assert readType in ['single', 'paired', 'combined']
	assert genome in ['host', 'virus']

	merge = bool(get_value_from_df(wildcards, 'merge'))
	dedup = bool(get_value_from_df(wildcards, 'dedup'))
	
	# if we want host alignment
	if genome == "virus":
		# if we're doing deduplication
		if dedup is True:
			# if we want single reads
			if readType == "single":
				return path.splitext(rules.align_bwa_virus.output.single)[0] + ".rmdup.sam"
			# if we want paired reads
			elif readType == 'paired':
				return path.splitext(rules.align_bwa_virus.output.paired)[0] + ".rmdup.sam"
			# if we want combined reads
			else:
				return path.splitext(rules.align_bwa_virus.output.combined)[0] + ".rmdup.sam"
		# if we're not doing deduplication
		else:
			if readType == "single":
				return rules.align_bwa_virus.output.single
			elif readType == "paired":
				return rules.align_bwa_virus.output.paired
			else:
				return rules.align_bwa_virus.output.combined
	# if we want the host alignment
	else:
		# if we're doing deduplication
		if dedup is True:
			# if we want single reads
			if readType == "single":
				return path.splitext(rules.align_bwa_host_single.output.sam)[0] + ".rmdup.sam"
			elif readType == "paired":
				return path.splitext(rules.align_bwa_host_paired.output.sam)[0] + ".rmdup.sam"
			else:
				return path.splitext(rules.combine_host.output.combined)[0] + ".rmdup.sam"
		# if we're not doing deduplication
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
		fastq = temp("{outpath}/{dset}/virus_aligned/{samp}.bwaSingle.mappedTo{virus}.fastq.gz"),
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
		samtools view -h -F 0x4 -F 0x800 -F 0x100 -o - {input.aligned} | samtools fastq -0 {output.fastq} - 
		"""

rule extract_vAligned_paired:
	input:
		aligned = lambda wildcards: get_sam(wildcards, "paired", "virus")
	output:
		pvBam_readMap_mateUnmap = temp("{outpath}/{dset}/virus_aligned/{samp}.{virus}.bwaPaired.L.bam"),
		pvBam_readUnmap_mateMap = temp("{outpath}/{dset}/virus_aligned/{samp}.{virus}.bwaPaired.R.bam"),
		pvBam_bothMapped = temp("{outpath}/{dset}/virus_aligned/{samp}.{virus}.bwaPaired.B.bam"),
		bam = temp("{outpath}/{dset}/virus_aligned/{samp}.{virus}.bwaPaired.mapped.bam"),
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
		fastq1 = temp("{outpath}/{dset}/virus_aligned/{samp}.bwaPaired.mappedTo{virus}.1.fastq.gz"),
		fastq2 = temp("{outpath}/{dset}/virus_aligned/{samp}.bwaPaired.mappedTo{virus}.2.fastq.gz")
	group: "host_align"
	conda:
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	shell:
		"""
		samtools collate -O {input.bam} | samtools fastq -1 {output.fastq1} -2 {output.fastq2} -
		"""

rule align_bwa_host_single:
	input:	
		idx = expand("{outpath}/references/{host}/{host}.{ext}", ext=["ann", "amb", "bwt", "pac", "sa"], allow_missing=True),
		all = rules.extract_to_fastq_single.output.fastq,
	output:
		sam = temp("{outpath}/{dset}/host_aligned/{samp}.{host}.readsFrom{virus}.bwaSingle.sam"),
	group: "host_align"
	conda: 
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	resources:
		mem_mb=lambda wildcards, attempt, input: attempt * 3 * sum([int(os.stat(file).st_size/1e6) for file in input.idx])
	params:
		index = lambda wildcards, input: path.splitext(input.idx[0])[0],
		mapping = lambda wildcards: get_value_from_df(wildcards, 'bwa_mem_params'),
		RG = lambda wildcards: f"-R '@RG\\tID:{wildcards.samp}_{wildcards.host}_merged\\tSM:{wildcards.samp}\\tPM:merged'"
	threads: 4
	shell:		
		"""
		bwa mem -t {threads} {params.mapping} {params.RG} -o {output.sam} {params.index} {input.all}
		"""
		
rule align_bwa_host_paired:
	input:	
		idx = expand("{outpath}/references/{host}/{host}.{ext}", ext=["ann", "amb", "bwt", "pac", "sa"], allow_missing=True),
		r1 = rules.extract_to_fastq_paired.output[0],
		r2 = rules.extract_to_fastq_paired.output[1]
	output:
		sam = temp("{outpath}/{dset}/host_aligned/{samp}.{host}.readsFrom{virus}.bwaPaired.sam"),
	group: "host_align"
	conda: 
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	resources:
		mem_mb=lambda wildcards, attempt, input: attempt * 3 * sum([int(os.stat(file).st_size/1e6) for file in input.idx])
	params:
		index = lambda wildcards, input: path.splitext(input.idx[0])[0],
		mapping = lambda wildcards: get_value_from_df(wildcards, 'bwa_mem_params'),
		RG = lambda wildcards: f"-R '@RG\\tID:{wildcards.samp}_{wildcards.host}_unmerged\\tSM:{wildcards.samp}\\tPM:unmerged'"
	threads: 4
	shell:		
		"""
		bwa mem -t {threads} {params.mapping} {params.RG} -o {output.sam} {params.index} {input.r1} {input.r2}
		"""
		
rule combine_host:
	input:
		paired = rules.align_bwa_host_paired.output.sam,
		single = rules.align_bwa_host_single.output.sam
	output:
		combined = temp("{outpath}/{dset}/host_aligned/{samp}.{host}.readsFrom{virus}.sam")
	group: "host_align"
	conda: 
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	resources:
		mem_mb=lambda wildcards, attempt, input: attempt * 3 * sum([int(os.stat(file).st_size/1e6) for file in input])
	shell:		
		"""
		samtools merge {output.combined} {input.single} {input.paired}
		"""
	

#### sam file manipulations ####

rule convert_to_bam:
	input:
		sam = "{outpath}/{dset}/{folder}/{alignment}.sam"
	output:
		bam = "{outpath}/{dset}/{folder}/{alignment}.bam",
		bai = "{outpath}/{dset}/{folder}/{alignment}.bam.bai"
	params:
		tmp_prefix = lambda wildcards, input: path.splitext(input[0])[0]
	wildcard_constraints:
		folder = "host_aligned|virus_aligned"
	conda: 
		"../envs/bwa.yml"	
	container:
		"docker://szsctt/bwa:1"
	resources:
		mem_mb=lambda wildcards, attempt, input: attempt * 3 * int(os.stat(input.sam).st_size/1e6),
	shell:
		"""
		rm -f {params.tmp_prefix}*tmp*
		samtools sort -o {output.bam} {input.sam}
		samtools index {output.bam}
		"""
		
rule markdup:
	input:
		sam = "{outpath}/{dset}/{folder}/{alignment}.sam"
	output:
		fixmate = temp("{outpath}/{dset}/{folder}/{alignment}.fixmate.bam"),
		markdup = temp("{outpath}/{dset}/{folder}/{alignment}.dups.sam"),
		metrics = temp("{outpath}/{dset}/{folder}/{alignment}.dups.txt")
	wildcard_constraints:
		folder = "host_aligned|virus_aligned"
	group: "rmdup"
	resources:
		mem_mb=lambda wildcards, attempt, input: attempt * 3 * int(os.stat(input.sam).st_size/1e6)
	conda: 
		"../envs/picard.yml"	
	container:
		"docker://szsctt/picard:1"
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
	group: "rmdup"
	wildcard_constraints:
		folder = "host_aligned|virus_aligned"
	conda: 
		"../envs/bwa.yml"	
	container:
		"docker://szsctt/bwa:1"
	resources:
		mem_mb=lambda wildcards, attempt, input: attempt * 3 * int(os.stat(input.sam).st_size/1e6)
	shell:
		"""
		samtools view -h -F 1024 {input.sam} > {output.sam}
		"""
	

