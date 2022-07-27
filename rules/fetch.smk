def get_acc_list(wildcards):
	dir = checkpoints.load_tsv.get(**wildcards).output[0]
	accs = glob_wildcards(os.path.join(dir, "{acc}.txt")).acc
	if len(accs) == 0:
		raise ValueError("No accessions found!")
	return accs

		
rule fq_dump:
	input:
		tsv = rules.efetch.output.xml
	output:
		r1 = "out/sra_fastq/{acc}_1.fastq",
		r2 = "out/sra_fastq/{acc}_2.fastq",
		unpaired = "out/sra_fastq/{acc}.fastq"
	conda:
		"envs/sra_toolkit.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/fq_dump/{acc}.log"
	params:
		outdir = lambda wildcards, output: f"-O {os.path.dirname(output.r1)}",
	resources:
		disk_mb = lambda wildcards, attempt: 10000 + (10000 * attempt),
		mem_mb = lambda wildcards, attempt: 2000 * attempt
	threads: 4
	shell:
		"""
		fasterq-dump --temp {resources.tmpdir} --threads {threads} \
		--mem {resources.mem_mb} {params} {wildcards.acc}
		touch {output}
		"""

rule bbduk_single:
	input:
		unpaired = rules.fq_dump.output.unpaired,
		adapters = "data/adapters.fa"
	output:
		unpaired = "out/bbduk_single/{acc}.fastq"
	conda:
		"envs/bbtools.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/bbduk_single/{acc}.log"	
	params:
		eoom = "-eoom",
		mem = lambda wildcards, resources: f"-Xmx{int(resources.mem_mb)}m"	
	resources:
		disk_mb = lambda wildcards, attempt: 10000 + (10000 * attempt),
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb, 2000) * 5
	shell:
		"""
		bbduk.sh {params} in={input} out={output} ref={input.adapters} ktrim=r k=23 mink=11 hdist=1
		"""

rule bbduk_paired:
	input:
		r1 = rules.fq_dump.output.r1,
		r2 = rules.fq_dump.output.r2,
		adapters = "data/adapters.fa"
	output:
		r1 = "out/bbduk_paired/{acc}_1.fastq",
		r2 = "out/bbduk_paired/{acc}_2.fastq",
	conda:
		"envs/bbtools.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/bbduk_paired/{acc}.log"	
	params:
		eoom = "-eoom",
		mem = lambda wildcards, resources: f"-Xmx{int(resources.mem_mb)}m"	
	resources:
		disk_mb = lambda wildcards, attempt: 10000 + (10000 * attempt),
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb * 5, 2000) * attempt
	shell:
		"""
		bbduk.sh {params} in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} \
		ref={input.adapters} ktrim=r k=23 mink=11 hdist=1 tpe tbo
		"""
