rule index_caps:
	input:
		fa = small_nt_refs
	output:
		 idx=multiext("out/bwa_idx/caps", ".amb", ".ann", ".bwt", ".pac", ".sa"),
	conda:
		"envs/bwa.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/bwa_index/bwa_index.log"
	params:
		prefix = lambda wildcards, output: f"-p {os.path.splitext(output.idx[0])[0]}"
	resources:
		disk_mb = lambda wildcards, attempt, input: max(input.size_mb * 5, 5000) * attempt,
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb * 5, 2000)	* attempt
	shell:
		"""
		bwa index {params} {input.fa}
		"""

rule remap_bwa_single:
	input:
		unpaired = rules.extract_single_reads.output.fastq,
		idx = rules.index_caps.output.idx
	output:
		bam = "out/remap_bwa_single/{acc}.bam"
	conda:
		"envs/bwa.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/remap_bwa_single/{acc}.log"
	params:
		prefix = lambda wildcards, input: os.path.splitext(input.idx[0])[0],
		score_output = '-T 10',
		read_group = lambda wildcards: f"-R '@RG\\tID:{wildcards.acc}_single_remap\\tSM:{wildcards.acc}'"
	resources:
		disk_mb = lambda wildcards, attempt, input: max(input.size_mb * 5, 10000) * attempt,
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb * 3, 2000)	* attempt
	threads: 8
	shell:
		"""
		bwa mem -t {threads} {params.score_output} {params.read_group}\
            {params.prefix} {input.unpaired} |\
		samtools sort -o {output.unpaired} - 
		samtools index {output.unpaired}
		"""	
		
rule remap_bwa_paired:
	input:
		r1 = rules.bbduk_paired.output.r1,
		r2 = rules.bbduk_paired.output.r2,
		idx = rules.index_caps.output.idx
	output:
		aligned = "out/remap_bwa_paired/{acc}.bam",
		idx = "out/remap_bwa_paired/{acc}.bam.bai"
	conda:
		"envs/bwa.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/remap_bwa_paired/{acc}.log"
	params:
		prefix = lambda wildcards, input: os.path.splitext(input.idx[0])[0],
		score_output = '-T 10',
		read_group = lambda wildcards: f"-R '@RG\\tID:{wildcards.acc}_paired_remap\\tSM:{wildcards.acc}'"
	resources:
		disk_mb = lambda wildcards, attempt, input: max(input.size_mb * 5, 10000) * attempt,
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb * 3, 2000)	* attempt
	threads: 8
	shell:
		"""
		bwa mem -t {threads} {params.score_output} {params.read_group} \
            {params.prefix} {input.r1} {input.r2} |\
		samtools sort -o {output.aligned} -
		samtools index {output.aligned}
		"""	
		
		
rule merge_remap:
	input:
		single = rules.remap_bwa_single,
		paried = rules.remap_bwa_paired,
	output:
		merged = "out/remap_bwa/{acc}.bam"
	conda:
		"envs/bwa.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logsremap_bwa_paired/{acc}.log"
	resources:
		disk_mb = lambda wildcards, attempt, input: max(input.size_mb * 5, 10000) * attempt,
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb * 3, 2000)	* attempt
	shell:
		"""
		samtools merge {output} {input}
		"""		
		
rule call_cap_variants:
	input:
		merged = rules.merge_remap.output.merged,
		ref = rules.index_caps.input.fa
	output:
		vcf = "out/freebayes/{acc}.vcf"
	conda:
		"envs/freebayes.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logsremap_bwa_paired/{acc}.log"
	resources:
		disk_mb = lambda wildcards, attempt, input: max(input.size_mb * 5, 10000) * attempt,
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb * 3, 2000)	* attempt
	shell:
		"""
		freebayes -f {input.ref} {input.merged} > {output.vcf}
		"""