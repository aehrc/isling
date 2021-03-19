
rule run_soft:
	input:
		host = lambda wildcards: get_sam(wildcards, "combined", "host"),
		virus = lambda wildcards: get_sam(wildcards, "combined", "virus")
	output:
		soft = temp("{outpath}/{dset}/ints/{samp}.{part}.{host}.{virus}.soft.txt"),
	params:
		cutoff = lambda wildcards: f"--cutoff {int(get_value_from_df(wildcards, 'clip_cutoff'))}",
		tol = lambda wildcards: f"--tol {int(get_value_from_df(wildcards, 'cigar_tol'))}",
	resources:
		mem_mb=lambda wildcards, attempt, input: int(resources_list_with_min_and_max(input, attempt, 1.5)),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
		nodes = 1
	container:
		"docker://ubuntu:18.04"	
	group: "ints"
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
	params:
		cutoff = lambda wildcards: f"--cutoff {int(get_value_from_df(wildcards, 'clip_cutoff'))}",
		tol = lambda wildcards: f"--tol {int(get_value_from_df(wildcards, 'cigar_tol'))}",
	resources:
		mem_mb=lambda wildcards, attempt, input: int(resources_list_with_min_and_max(input, attempt, 1.5)),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
	container:
		"docker://ubuntu:18.04"
	group: "ints"
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
	params:
		cutoff = lambda wildcards: f"--cutoff {int(get_value_from_df(wildcards, 'clip_cutoff'))}",
		tol = lambda wildcards: f"--tol {int(get_value_from_df(wildcards, 'cigar_tol'))}",
		tlen = lambda wildcards: f"--tlen {get_value_from_df(wildcards, 'mean_frag_len')}"
	resources:
		mem_mb=lambda wildcards, attempt, input: int(resources_list_with_min_and_max(input, attempt, 1.5)),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
	container:
		"docker://ubuntu:18.04"
#	threads: workflow.cores
	group: "ints"
	shell:
		"""
		perl -Iscripts scripts/discordant.pl --viral {input.virus} --host {input.host} --output {output.discord} {params}
		"""

rule combine_ints:
	input:
		soft = lambda wildcards: expand(strip_wildcard_constraints(rules.run_soft.output.soft), 
										part = get_split(wildcards), allow_missing = True),
		short = lambda wildcards: expand(strip_wildcard_constraints(rules.run_short.output.short), 
										part = get_split(wildcards), allow_missing = True),
		discordant = lambda wildcards: expand(strip_wildcard_constraints(rules.run_discordant.output.discord), 
										part = get_split(wildcards), allow_missing = True)
	output:
		all = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.txt"
	resources:
		mem_mb=lambda wildcards, attempt, input: int(resources_list_with_min_and_max(input, attempt, 1.5)),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
	container:
		"docker://ubuntu:18.04"
	group: "ints"
	shell:
		"""
		(head -n1 {input.soft[0]} && awk '(FNR==1){{next}}{{print $0| "sort -k1,1 -k2,2n"}}' {input}) > {output.all}
		"""

