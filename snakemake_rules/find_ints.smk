
rule run_soft:
	input:
		host = lambda wildcards: get_sam(wildcards, "combined", "host"),
		virus = lambda wildcards: get_sam(wildcards, "combined", "virus")
	output:
		soft = temp("{outpath}/{dset}/ints/{samp}.{part}.{host}.{virus}.soft.txt"),
	group: "ints"
	params:
		cutoff = lambda wildcards: f"--cutoff {int(get_value_from_df(wildcards, 'clip_cutoff'))}",
		tol = lambda wildcards: f"--tol {int(get_value_from_df(wildcards, 'cigar_tol'))}",
		min_mapq = lambda wildcards: f"--min-mapq {int(get_value_from_df(wildcards, 'min_mapq'))}",
	resources:
		mem_mb=lambda wildcards, attempt, input: resources_list_with_min_and_max((input.host, input.virus), attempt)
	container:
		"docker://ubuntu:18.04"	
	shell:
		"""
		perl -Iscripts scripts/softClip.pl --viral {input.virus} --host {input.host} --output {output.soft} {params}
		"""
		
rule run_short:
	input:
		host = lambda wildcards: get_sam(wildcards, "combined", "host"),
		virus = lambda wildcards: get_sam(wildcards, "combined", "virus"),
	output:
		short = temp("{outpath}/{dset}/ints/{samp}.{part}.{host}.{virus}.short.txt"),
	group: "ints"
	params:
		cutoff = lambda wildcards: f"--cutoff {int(get_value_from_df(wildcards, 'clip_cutoff'))}",
		tol = lambda wildcards: f"--tol {int(get_value_from_df(wildcards, 'cigar_tol'))}",
		min_mapq = lambda wildcards: f"--min-mapq {int(get_value_from_df(wildcards, 'min_mapq'))}",
	resources:
		mem_mb=lambda wildcards, attempt, input: resources_list_with_min_and_max((input.host, input.virus), attempt)
	container:
		"docker://ubuntu:18.04"
	shell:
		"""
		perl -Iscripts scripts/short.pl --viral {input.virus} --host {input.host} --output {output.short} {params}
		"""
		
rule run_discordant:
	input:
		host = lambda wildcards: get_sam(wildcards, "paired", "host"),
		virus = lambda wildcards: get_sam(wildcards, "paired", "virus"),
	output:
		discord = temp("{outpath}/{dset}/ints/{samp}.{part}.{host}.{virus}.discordant.txt"),
	group: "ints"
	params:
		cutoff = lambda wildcards: f"--cutoff {int(get_value_from_df(wildcards, 'clip_cutoff'))}",
		tol = lambda wildcards: f"--tol {int(get_value_from_df(wildcards, 'cigar_tol'))}",
		min_mapq = lambda wildcards: f"--min-mapq {int(get_value_from_df(wildcards, 'min_mapq'))}",
		tlen = lambda wildcards: f"--tlen {get_value_from_df(wildcards, 'mean_frag_len')}"
	resources:
		mem_mb=lambda wildcards, attempt, input: resources_list_with_min_and_max((input.host, input.virus), attempt)
	container:
		"docker://ubuntu:18.04"
	threads: workflow.cores
	shell:
		"""
		perl -Iscripts scripts/discordant.pl --viral {input.virus} --host {input.host} --output {output.discord} {params}
		"""

rule combine_ints:
	input:
		soft = rules.run_soft.output,
		short = rules.run_short.output,
		discordant = rules.run_discordant.output
	group: "ints"
	output:
		all = temp("{outpath}/{dset}/ints/{samp}.{part}.{host}.{virus}.integrations.txt")
	container:
		"docker://ubuntu:18.04"
	shell:
		"""
		awk 'FNR>1 || NR==1' {input} > {output.all}
		"""

rule merge_parts_ints:
	message: "Merge the intergrations from each part into one file"
	input:
		files = expand("{{outpath}}/{{dset}}/ints/{{samp}}.{parts}.{{host}}.{{virus}}.integrations.txt", parts = get_split()),
	group: "ints"
	output:
		all = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.txt",
		temp =  temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.txt.tmp"),
	container:
		"docker://ubuntu:18.04"
	threads: 1
	shell:
		"""
		echo {input.files}
		awk 'FNR>1 || NR==1' {input.files} > {output.all}
		sort -n -k1,1 -k2,2n {output.all} > {output.temp}
		cp {output.temp} {output.all}
		"""
