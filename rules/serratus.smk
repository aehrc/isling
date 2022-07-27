
rule fetch_serratus:
	output:
		aav = "out/serratus/AAV_test.tsv",
		parv = "out/serratus/parvoviridae.tsv"
	conda:
		"envs/R.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/fetch_serratus/log"	
	shell:
		"""
		Rscript scripts/get_sra_accs_from_serratus.R -aav {output.aav} -parv {output.parv}
		"""	

checkpoint load_tsv:
	input:
		tsv = rules.fetch_serratus.output.aav
	output:
		directory("out/accs/")
	conda:
		"envs/parallel.yml"
	shell:
		"""
		ACCS=$(awk -F'\t' 'NR>1 {{print $1}} ' {input.tsv})
		echo $ACCS
		mkdir -p {output}
		cd {output}
		parallel touch {{}}.txt ::: ${{ACCS}}
		"""

rule efetch:
	input:
		acc = "out/accs/{acc}.txt"
	output:
		xml = "out/info_xml/{acc}.xml"
	conda:
		"envs/eutils.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/efetch/{acc}.log"
	shell:
		"""
		esearch -db sra -query {wildcards.acc} | efetch -format xml > {output.xml}
		"""
		
rule parse_xml:
	input:
		xml = rules.efetch.output.xml
	output:
		tsv = "out/info_tsv/{acc}.tsv"
	conda:
		"envs/R.yml"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/parse_xml/{acc}.log"
	shell:
		"""
		time Rscript scripts/extract_xml_info.R \
			--input {input.xml} \
			--output {output.tsv} \
			--accession {wildcards.acc}
		"""

rule combine_tsvs:
	input:
		tsvs = lambda wildcards: expand("out/info_tsv/{acc}.tsv", acc=get_acc_list(wildcards))
	output:
		all = "out/info.tsv"
	singularity:
		"docker://szsctt/mine:latest"
	log:
		"logs/combine_info.log"
	shell:
		"""
		(head -n1 {input.tsvs[0]} && awk '(FNR==1){{next}}{{print}}' {input}) > {output.all}
		"""