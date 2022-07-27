
rule megahit:
	input:
		unpaired = rules.extract_single_reads.output.fastq,
		r1 = rules.extract_paired_reads.output.fq1,
		r2 = rules.extract_paired_reads.output.fq2
	output:
		contigs = "out/megahit/{acc}/{acc}.contigs.fa"
	conda:
		"envs/megahit.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/megahit/{acc}.log"
	params:
		dir = lambda wildcards, output: os.path.dirname(output.contigs),
		prefix = lambda wildcards: f"--out-prefix {wildcards.acc}",
		mem = lambda wildcards, resources: f"-m {resources.mem_mb * 1024 * 1024}"
	resources:
		disk_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 10000) * attempt,
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb * 6, 2000) * attempt
	threads: 8
	shell:
		"""
		rm -r {params.dir}
		megahit -1 {input.r1} -2 {input.r2} -r {input.unpaired} \
		-o {params} -t {threads} 
		"""
		
rule minimap:
	input:
		contigs = rules.megahit.output.contigs,
		ref = small_nt_refs
	output:
		bam = "out/minimap/{acc}.bam",
		idx = "out/minimap/{acc}.bam.bai"
	conda:
		"envs/minimap.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/minimap/{acc}.log"	
	resources:
		disk_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 5000) * attempt,
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 2000) * attempt
	threads: 4
	shell:
		"""
		minimap2 -ax asm20 -t {threads} {input.ref} {input.contigs} |\
		samtools sort -o {output.bam} -
		samtools index {output.bam}
		"""

rule minimap_flagstat:
	input:
		bam = rules.minimap.output.bam
	output:
		txt = "out/minimap_txt/{acc}.txt"
	conda:
		"envs/minimap.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/minimap_flagstat/{acc}.log"	
	resources:
		disk_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 5000) * attempt,
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 2000) * attempt
	shell:
		"""
		samtools flagstat {input.bam} > {output.txt}
		"""	

rule done_minimap:
	input:
		lambda wildcards: expand("out/minimap_txt/{acc}.txt", acc=get_acc_list(wildcards))
	output:
		"done_minimap"
	shell:
		"""
		touch {output}
		"""	