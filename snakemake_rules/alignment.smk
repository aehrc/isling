		
#functions for if we did seqPrep or not
def get_for_align(wildcards, read_type):
	row_idx = list(toDo.loc[:,'unique']).index(f"{wildcards.dset}+++{wildcards.samp}")
	merge_check = bool(toDo.loc[toDo.index[row_idx], 'merge'])
	if merge_check is True:
		folder = "merged_reads"
	else:
		folder = "combined_reads"
	return f"{wildcards.outpath}/{wildcards.dset}/{folder}/{wildcards.samp}.{read_type}.fastq.gz"

def get_value_from_df(wildcards, column_name):
	# get a value from the row of the df corresponding to this sample and dataset
	row_idx = list(toDo.loc[:,'unique']).index(f"{wildcards.dset}+++{wildcards.samp}")
	return toDo.loc[toDo.index[row_idx], column_name]


#### alignments ####

rule index:
	input:
		lambda wildcards: ref_names[wildcards.genome]
	output:
		expand("{outpath}/references/{genome}/{genome}.{ext}", ext=["ann", "amb", "bwt", "pac", "sa"], allow_missing=True)
	conda: 
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	params:
		prefix = lambda wildcards, output: path.splitext(output[0])[0]
	shell:
		"bwa index -p {params.prefix} {input}"

	
rule align_bwa_virus_single:
	input:
		idx = expand("{outpath}/references/{virus}/{virus}.{ext}", ext=["ann", "amb", "bwt", "pac", "sa"], allow_missing=True),
		all = lambda wildcards: get_for_align(wildcards, "all"),
	output:
		sam = temp("{outpath}/{dset}/virus_aligned/{samp}.{virus}.bwaSingle.sam")
	params:
		index = lambda wildcards, input: path.splitext(input.idx[0])[0],
		mapping = lambda wildcards: get_value_from_df(wildcards, 'bwa_mem_params')
	resources:
		mem_mb=lambda wildcards, attempt: attempt * 4000
	conda:
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	threads: 5
	shell:
		"""
		bwa mem -t {threads} {params.mapping} -o {output.sam} {params.index} {input.all}
		"""

rule align_bwa_virus_paired:
	input:
		idx = expand("{outpath}/references/{virus}/{virus}.{ext}", ext=["ann", "amb", "bwt", "pac", "sa"], allow_missing=True),
		r1 = lambda wildcards: get_for_align(wildcards, "1"),
		r2 = lambda wildcards: get_for_align(wildcards, "2"),
	output:
		sam = temp("{outpath}/{dset}/virus_aligned/{samp}.{virus}.bwaPaired.sam")
	params:
		index = lambda wildcards, input: path.splitext(input.idx[0])[0],
		mapping = lambda wildcards: get_value_from_df(wildcards, 'bwa_mem_params')
	resources:
		mem_mb=lambda wildcards, attempt: attempt * 4000	
	conda:
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	threads: 5
	shell:
		"""
		bwa mem -t {threads} {params.mapping} -o {output.sam} {params.index} {input.r1} {input.r2} 
		"""
		

def get_sam(wildcards, readType, genome):
	row_idx = list(toDo.loc[:,'unique']).index(f"{wildcards.dset}+++{wildcards.samp}")
	dedup_check = bool(toDo.loc[toDo.index[row_idx], 'dedup'])
	# if we want host alignment
	if genome == "virus":
		# if we're doing deduplication
		if dedup_check is True:
			# if we want single reads
			if readType == "single":
				return path.splitext(rules.align_bwa_virus_single.output[0])[0] + ".rmdup.sam"
			else:
				return path.splitext(rules.align_bwa_virus_paired.output[0])[0] + ".rmdup.sam"
		# if we're not doing deduplication
		else:
			if readType == "single":
				return path.splitext(rules.align_bwa_virus_single.output[0])[0] + ".dups.sam"
			else:
				return path.splitext(rules.align_bwa_virus_paired.output[0])[0] + ".dups.sam"
	# if we want the host alignment
	else:
		# if we're doing deduplication
		if dedup_check is True:
			# if we want single reads
			if readType == "single":
				return path.splitext(rules.align_bwa_host_single.output[0])[0] + ".rmdup.sam"
			else:
				return path.splitext(rules.align_bwa_host_paired.output[0])[0] + ".rmdup.sam"
		# if we're not doing deduplication
		else:
			if readType == "single":
				return path.splitext(rules.align_bwa_host_single.output[0])[0] + ".dups.sam"
			else:
				return path.splitext(rules.align_bwa_host_paired.output[0])[0] + ".dups.sam"	
				
		
		
rule extract_vAligned_single:
	input:
		aligned = lambda wildcards: get_sam(wildcards, "single", "virus"),
	output:
		sam = temp("{outpath}/{dset}/virus_aligned/{samp}.{virus}.bwaSingle.mapped.sam"),
	conda:
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	shell:
		"""
		samtools view -h -F 0x4 -F 0x800 -o - {input.aligned} | samtools sort - -n -o {output.sam}
		"""

rule extract_vAligned_paired:
	input:
		aligned = lambda wildcards: get_sam(wildcards, "paired", "virus")
	output:
		pvBam_readMap_mateUnmap = temp("{outpath}/{dset}/virus_aligned/{samp}.{virus}.bwaPaired.L.bam"),
		pvBam_readUnmap_mateMap = temp("{outpath}/{dset}/virus_aligned/{samp}.{virus}.bwaPaired.R.bam"),
		pvBam_bothMapped = temp("{outpath}/{dset}/virus_aligned/{samp}.{virus}.bwaPaired.B.bam"),
		sam = temp("{outpath}/{dset}/virus_aligned/{samp}.{virus}.bwaPaired.mapped.sam"),
	conda:
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	shell:
		"""
		samtools view -hb -F 0x4 -f 0x8 -F 0x800 -o {output.pvBam_readMap_mateUnmap} {input.aligned}
		samtools view -hb -f 0x4 -F 0x8 -F 0x800 -o {output.pvBam_readUnmap_mateMap} {input.aligned}
		samtools view -hb -F 0x4 -F 0x8 -F 0x800 -o {output.pvBam_bothMapped} {input.aligned}
		samtools merge - {output.pvBam_readMap_mateUnmap} {output.pvBam_bothMapped} {output.pvBam_readUnmap_mateMap} | samtools sort - -n -o {output.sam}
		"""
		

rule extract_to_fastq_single:
	input:
		sam = rules.extract_vAligned_single.output.sam
	output:
		fastq = temp("{outpath}/{dset}/virus_aligned/{samp}.bwaSingle.mappedTo{virus}.fastq.gz")
	conda:
		"../envs/picard.yml"
	container:
		"docker://szsctt/picard:1"
	shell:
		"""
		picard SamToFastq I={input.sam} FASTQ={output.fastq}
		"""
		
rule extract_to_fastq_paired:
	input:
		sam = rules.extract_vAligned_paired.output.sam
	output:
		fastq1 = temp("{outpath}/{dset}/virus_aligned/{samp}.bwaSingle.mappedTo{virus}.1.fastq.gz"),
		fastq2 = temp("{outpath}/{dset}/virus_aligned/{samp}.bwaSingle.mappedTo{virus}.2.fastq.gz")
	conda:
		"../envs/picard.yml"
	container:
		"docker://szsctt/picard:1"
	shell:
		"""
		picard SamToFastq I={input.sam} FASTQ={output.fastq1} SECOND_END_FASTQ={output.fastq2}
		"""

rule align_bwa_host_single:
	input:	
		idx = expand("{outpath}/references/{host}/{host}.{ext}", ext=["ann", "amb", "bwt", "pac", "sa"], allow_missing=True),
		all = rules.extract_to_fastq_single.output.fastq
	output:
		sam = temp("{outpath}/{dset}/host_aligned/{samp}.{host}.readsFrom{virus}.bwaSingle.sam"),
	conda: 
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	resources:
		mem_mb=lambda wildcards, attempt: attempt * 4000
	params:
		index = lambda wildcards, input: path.splitext(input.idx[0])[0],
		mapping = lambda wildcards: get_value_from_df(wildcards, 'bwa_mem_params')
	threads: 4
	shell:		
		"""
		bwa mem -t {threads} {params.mapping} -o {output.sam} {params.index} {input.all}
		"""
		
rule align_bwa_host_paired:
	input:	
		idx = expand("{outpath}/references/{host}/{host}.{ext}", ext=["ann", "amb", "bwt", "pac", "sa"], allow_missing=True),
		r1 = rules.extract_to_fastq_paired.output[0],
		r2 = rules.extract_to_fastq_paired.output[1]
	output:
		sam = temp("{outpath}/{dset}/host_aligned/{samp}.{host}.readsFrom{virus}.bwaPaired.sam"),
	conda: 
		"../envs/bwa.yml"
	container:
		"docker://szsctt/bwa:1"
	resources:
		mem_mb=lambda wildcards, attempt: attempt * 4000	
	params:
		index = lambda wildcards, input: path.splitext(input.idx[0])[0],
		mapping = lambda wildcards: get_value_from_df(wildcards, 'bwa_mem_params')
	threads: 4
	shell:		
		"""
		bwa mem -t {threads} {params.mapping} -o {output.sam} {params.index} {input.r1} {input.r2}
		"""

#### sam file manipulations ####

rule convert_to_bam:
	input:
		"{outpath}/{dset}/{folder}/{alignment}.sam"
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
		mem_mb=lambda wildcards, attempt: attempt * 4000
	shell:
		"""
		rm -f {params.tmp_prefix}*tmp*
		samtools sort -o {output.bam} {input}
		samtools index {output.bam}
		"""
		
rule markdup:
	input:
		sam = "{outpath}/{dset}/{folder}/{alignment}.sam"
	output:
		fixmate = temp("{outpath}/{dset}/{folder}/{alignment}.fixmate.bam"),
		markdup = "{outpath}/{dset}/{folder}/{alignment}.dups.sam",
		metrics = temp("{outpath}/{dset}/{folder}/{alignment}.dups.txt")
	wildcard_constraints:
		folder = "host_aligned|virus_aligned"
	resources:
		mem_mb=lambda wildcards, attempt: attempt * 4000
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
	wildcard_constraints:
		folder = "host_aligned|virus_aligned"
	conda: 
		"../envs/bwa.yml"	
	container:
		"docker://szsctt/bwa:1"
	shell:
		"""
		samtools view -h -F 1024 {input.sam} > {output.sam}
		"""
	

