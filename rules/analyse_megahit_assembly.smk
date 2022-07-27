		
		
rule prokka:
	input:
		contigs = rules.megahit.output.contigs,
	output:
		prokka = multiext("out/prokka/{acc}/{acc}", ".log", ".err", ".txt", ".tsv", 
							".tbl", ".fsa", ".sqn", ".ffn", ".faa", ".fna", ".gbk", ".gff")
	conda:
		"envs/prokka.yml"
	singularity:
		"docker://staphb/prokka:latest"
	log:
		"logs/prokka/{acc}.log"	
	params:
		outdir = lambda wildcards, output: os.path.dirname(output.prokka[0])
	resources:
		disk_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 5000) * attempt,
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 2000) * attempt
	threads: 4
	shell:
		"""
		prokka --kingdom Virus --outdir {params.outdir} --force \
		--prefix {wildcards.acc}  --cpu {threads} {input.contigs} 
		"""	
		
rule diamond_db:
	input:
		ref = small_aa_refs
	output:
		ref = "out/diamond_db/caps.dmnd"
	conda:
		"envs/diamond.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/diamond_db/log.log"	
	params:
		prefix = lambda wildcards, output: f"-d {os.path.splitext(output.ref)[0]}"
	resources:
		disk_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 5000) * attempt,
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 2000) * attempt
	shell:
		"""
		diamond makedb --in {input.ref} {params}
		"""		
		
rule diamond_aln:
	input:
		faa = "out/prokka/{acc}/{acc}.faa",
		ref = rules.diamond_db.output.ref
	output:
		aln = "out/diamond/{acc}.tsv"
	conda:
		"envs/diamond.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/diamond_aln/{acc}.log"	
	params:
		prefix = lambda wildcards, input: f"-d {os.path.splitext(input.ref)[0]}",
		masking = "--masking 0",
		outfmt = "--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart"
				" qend sstart send evalue bitscore nident mismatch gapopen gaps btop cigar"
	resources:
		disk_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 5000) * attempt,
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 2000) * attempt
	shell:
		"""
		diamond blastp {params} -q {input.faa} -o {output.aln}
		"""		

rule check_diamond_aln:
	input:
		aln = rules.diamond_aln.output.aln,
		caps = rules.diamond_db.input.ref
	output:
		fl = "out/diamond_fl/{acc}.tsv"
	conda:
		"envs/biopython.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/check_diamond_aln/{acc}.log"	
	resources:
		disk_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 5000) * attempt,
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 2000) * attempt
	shell:
		"""
		python3 scripts/check_diamond_for_fl.py \
		 --caps-fa {input.caps} --diamond-output {input.aln} --outfile {output.fl} \
		 --sra-run {wildcards.acc}
		"""		

rule combine_diamond_fl:
	input:
		lambda wildcards: expand("out/diamond_fl/{acc}.tsv", acc=get_acc_list(wildcards))
	output:
		tsv = "out/diamond_fl.tsv"
	log:
		"logs/combine_diamond_fl/log.log"		
	resources:
		disk_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 5000) * attempt,
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 2000) * attempt
	shell:
		"""
		cat 'NR != 1 || FNR == 1' {input} > {output}
		"""		

rule subset_prokka_fa:
	input:
		fl_tsv = rules.check_diamond_aln.output.fl,
		prokka_faa = rules.diamond_aln.input.faa
	output:
		fa = "out/prokka_fl/{acc}.fa"	
	conda:
		"envs/biopython.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/subset_prokka_fa/{acc}.log"	
	resources:
		disk_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 5000) * attempt,
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 2000) * attempt
	shell:
		"""
		python3 scripts/subset_fasta.py \
		 --prokka-fa {input.prokka_faa} --fl-tsv {input.fl_tsv} --outfile {output.fa}
		"""		
		
rule identify_variants:
	input:
		query = rules.prokka.output[7],
		ref = small_nt_refs,
		aln = rules.check_diamond_aln.output.fl
	output:
		variants = "out/variants/{acc}.tsv"
	conda:
		"envs/biopython.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/identify_variants/{acc}.log"	
	resources:
		disk_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 5000) * attempt,
		mem_mb = lambda wildcards, attempt, input: max(input.size_mb * 4, 2000) * attempt
	shell:
		"""
		python3 scripts/identify_variants.py\
		 --dmd {input.aln} --refs {input.ref} \
		 --queries {input.query} --output {output.variants}
		"""		
