
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
		nm_pc = lambda wildcards: "" if get_value_from_df(wildcards, 'nm_pc') is None else f"--nm-pc {get_value_from_df(wildcards, 'nm_pc')}",
		nm_diff = lambda wildcards: "" if get_value_from_df(wildcards, 'nm_diff') is None else f"--nm-diff {int(get_value_from_df(wildcards, 'nm_diff'))}",
		frag_len = lambda wildcards: get_value_from_df(wildcards, 'mean_frag_len')
	resources:
		mem_mb=lambda wildcards, attempt, input: int(resources_list_with_min_and_max(input, attempt, 1.5)),
		runtime = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
		nodes = 1
	container:
		"docker://szsctt/isling:latest"
	conda:
		"../envs/find_ints.yml"
	shell:
		"""
		if [ '{params.frag_len}' = 'estimate' ]; then
			HOST=$(perl -ne '/^SN\sinsert size average:\t(\d+\.\d+)/ && print "$1"' {input.host_stats})
			VIRUS=$(perl -ne '/^SN\sinsert size average:\t(\d+\.\d+)/ && print "$1"' {input.virus_stats})		
			FRAG=$(echo $HOST $VIRUS | awk '{{print ($1 + $2) / 2}}')
		else
			FRAG={params.frag_len}
		fi
		
		echo "using mean fragment length $FRAG"
		
		python3 scripts/find_ints.py \
		--host {input.host} \
		--virus {input.virus}\
		--mean-template-length $FRAG \
		--integrations {output.ints} {params.cutoff} {params.tol} {params.nm_pc} {params.nm_diff}
		"""

rule combine_ints:
	input:
		ints = lambda wildcards: expand(strip_wildcard_constraints(rules.find_ints.output.ints), 
										part = get_split(wildcards), allow_missing = True),
	output:
		all = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.txt"
	resources:
		mem_mb=lambda wildcards, attempt, input: int(resources_list_with_min_and_max(input, attempt, 1.5)),
		runtime = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
	container:
		"docker://szsctt/isling:latest"
	shell:
		"""
		(head -n1 {input.ints[0]} && awk '(FNR==1){{next}}{{print $0| "sort -k1,1 -k2,2n"}}' {input}) > {output.all}
		"""

