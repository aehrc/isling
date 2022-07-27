
rule bwa_index:
	input:
		aavs = large_nt_refs
	output:
		 idx=multiext("out/bwa_idx/aavs", ".amb", ".ann", ".bwt", ".pac", ".sa"),
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
		bwa index {params} {input.aavs}
		"""
		
rule bwa_single:
	input:
		unpaired = rules.bbduk_single.output.unpaired,
		idx = rules.bwa_index.output.idx
	output:
		unpaired = "out/bwa_single/{acc}.bam",
		idx = "out/bwa_single/{acc}.bam.bai"
	conda:
		"envs/bwa.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/bwa_single/{acc}.log"
	params:
		prefix = lambda wildcards, input: os.path.splitext(input.idx[0])[0],
		score_output = '-T 10',
		read_group = lambda wildcards: f"-R '@RG\\tID:{wildcards.acc}_single\\tSM:{wildcards.acc}'"
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
		
rule bwa_paired:
	input:
		r1 = rules.bbduk_paired.output.r1,
		r2 = rules.bbduk_paired.output.r2,
		idx = rules.bwa_index.output.idx
	output:
		aligned = "out/bwa_paired/{acc}.bam",
		idx = "out/bwa_paired/{acc}.bam.bai"
	conda:
		"envs/bwa.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/bwa_paired/{acc}.log"
	params:
		prefix = lambda wildcards, input: os.path.splitext(input.idx[0])[0],
		score_output = '-T 10',
		read_group = lambda wildcards: f"-R '@RG\\tID:{wildcards.acc}_paired\\tSM:{wildcards.acc}'"
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

rule flagstat_single:
	input:
		bam = rules.bwa_single.output.unpaired
	output:
		txt = "out/bwa_single_flagstat/{acc}.txt"
	conda:
		"envs/bwa.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/bwa_single_flagstat/{acc}.log"	
	resources:
		disk_mb = lambda wildcards, attempt, input: max(input.size_mb * 2, 5000) * attempt,
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb * 0.5, 2000) * attempt
	shell:
		"""
		samtools flagstat {input.bam} > {output.txt}
		"""		
		
rule flagstat_paired:
	input:
		bam = rules.bwa_paired.output.aligned
	output:
		txt = "out/bwa_paired_flagstat/{acc}.txt"
	conda:
		"envs/bwa.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/bwa_paired_flagstat/{acc}.log"	
	resources:
		disk_mb = lambda wildcards, attempt, input: max(input.size_mb * 2, 5000) * attempt,
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb * 0.5, 2000) * attempt
	shell:
		"""
		samtools flagstat {input.bam} > {output.txt}
		"""			
		
rule extract_single_reads:
	input:
		bam = rules.bwa_single.output.unpaired
	output:
		fastq = "out/bwa_single_fastq/{acc}.fastq.gz"
	conda:
		"envs/bwa.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/extract_single_reads/{acc}.log"	
	resources:
		disk_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 5000) * attempt,
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb * 1.5, 2000) * attempt
	shell:
		"""
		samtools view -h -F 0x4 -F 0x100 -F 0x800 {input.bam} | \
		samtools fastq -0 {output.fastq} -
		"""	
		
rule extract_paired_reads:
	input:
		bam = rules.bwa_paired.output.aligned
	output:
		fq1 = "out/bwa_paired_fastq/{acc}_R1.fastq.gz",
		fq2 = "out/bwa_paired_fastq/{acc}_R2.fastq.gz",
	conda:
		"envs/bwa.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/extract_paired_reads/{acc}.log"	
	resources:
		disk_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 5000) * attempt,
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb * 1.5, 2000) * attempt
	shell:
		"""
		# note that this will only extract reads where both mates are mapped
		samtools view -h -F 0x4 -F 0x8 -F 0x100 -F 0x800 {input.bam} | \
		 samtools collate -O - |\
		 samtools fastq -1 {output.fq1} -2 {output.fq2} -
		"""	
		