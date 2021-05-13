
rule find_ints:
	input:
		host = lambda wildcards: get_sam(wildcards, "combined", "host"),
		virus = lambda wildcards: get_sam(wildcards, "combined", "virus")
	output:
		ints = temp("{outpath}/{dset}/ints/{samp}.{part}.{host}.{virus}.txt"),
	params:
		cutoff = lambda wildcards: f"--map-thresh {int(get_value_from_df(wildcards, 'clip_cutoff'))}",
		tol = lambda wildcards: f"--tolerance {int(get_value_from_df(wildcards, 'cigar_tol'))}",
	resources:
		mem_mb=lambda wildcards, attempt, input: int(resources_list_with_min_and_max(input, attempt, 1.5)),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
		nodes = 1
	container:
		"docker://simvi:2"
	shell:
		"""
		python3 scripts/find_ints.py --host {input.host} --virus {input.virus} --integrations {output.ints} {params}
		"""

rule combine_ints:
	input:
		ints = lambda wildcards: expand(strip_wildcard_constraints(rules.find_ints.output.ints), 
										part = get_split(wildcards), allow_missing = True),
	output:
		all = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.txt"
	resources:
		mem_mb=lambda wildcards, attempt, input: int(resources_list_with_min_and_max(input, attempt, 1.5)),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
	container:
		"docker://ubuntu:18.04"
	shell:
		"""
		(head -n1 {input.ints[0]} && awk '(FNR==1){{next}}{{print $0| "sort -k1,1 -k2,2n"}}' {input}) > {output.all}
		"""

