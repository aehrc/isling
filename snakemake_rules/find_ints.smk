
rule find_ints:
	input:
		host = lambda wildcards: get_sam(wildcards, "combined", "host"),
		virus = lambda wildcards: get_sam(wildcards, "combined", "virus"),
		host_stats = rules.host_stats.output.stats,
		virus_stats =  rules.virus_stats.output.stats,
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
		"docker://szsctt/simvi:2"
	shell:
		"""
		HOST=$(perl -ne '/^SN\sinsert size average:\t(\d+\.\d+)/ && print "$1"' {input.host_stats})
		VIRUS=$(perl -ne '/^SN\sinsert size average:\t(\d+\.\d+)/ && print "$1"' {input.virus_stats})		
		FRAG=$(echo $HOST $VIRUS | awk '{{print ($1 + $2) / 2}}')
		python3 scripts/find_ints.py \
		--host {input.host} \
		--virus {input.virus}\
		--mean-template-length $FRAG \
		--integrations {output.ints} {params}
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

