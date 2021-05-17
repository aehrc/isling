rule post_filter:
	input:
		ints = rules.combine_ints.output.all,
	output:
		kept = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.filter.txt"),
		excluded = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.removed.txt"),
	params:
		filterstring = lambda wildcards: get_value_from_df(wildcards, 'filter')
	resources:
		mem_mb=lambda wildcards, attempt, input: int(resources_list_with_min_and_max(input, attempt, 1.5)),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
	container:
		"docker://szsctt/simvi:1"
	conda: "../envs/filter.yml"
	shell:
		"""
		python3 scripts/filter.py -i {input.ints} -k {output.kept} -e {output.excluded} -c '{params.filterstring}'
		"""

rule sort_bed:
	input:
		unsorted = "{name}.bed"
	output:
		sorted = "{name}.sorted.bed"
	container:
		"docker://ubuntu:18.04"
	shell:
		"sort -k1,1 -k2,2n {input.unsorted} > {output.sorted}"

rule exclude_bed:
	input:
		beds = lambda wildcards: [os.path.splitext(i)[0] + '.sorted.bed' for i in get_value_from_df(wildcards, 'bed_exclude')],
		filt = rules.post_filter.output.kept,
		excluded = rules.post_filter.output.excluded
	output:
		tmp = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.filter2.txt.tmp"),
		kept = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.filter2.txt"),
		excluded = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.removed2.txt"),
	resources:
		mem_mb=lambda wildcards, attempt, input: int(resources_list_with_min_and_max(input, attempt, 1.5)),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
	container:
		"docker://szsctt/bedtools:1"
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		cp {input.filt} {output.tmp}
		ARRAY=( {input.beds} )
		for bed in "${{ARRAY[@]}}"; do
			echo "excluding integrations that intersect with $bed"
			head -n1 {input.filt} > {output.kept}
			bedtools intersect -v -a {output.tmp} -b $bed >> {output.kept}
			bedtools intersect -u -a {output.tmp} -b $bed >> {output.excluded}
			cp {output.kept} {output.tmp}
		done
		"""

def get_for_include_bed(wildcards, file_type):
	assert file_type in {'kept', 'excluded'}
	# if there aren't any bed files to use for excluding
	if len(get_value_from_df(wildcards, 'bed_exclude')) == 0:
		if file_type == 'kept':
			return rules.post_filter.output.kept
		else:
			return rules.post_filter.output.excluded
	# if there are, then use files after excluding
	else:
		if file_type == 'kept':
			return rules.exclude_bed.output.kept
		else:
			return rules.exclude_bed.output.excluded
		

rule include_bed:
	input:
		beds = lambda wildcards: [os.path.splitext(i)[0] + '.sorted.bed' for i in get_value_from_df(wildcards, 'bed_include')],	
		filt = lambda wildcards:  get_for_include_bed(wildcards, 'kept'),
		excluded = lambda wildcards:  get_for_include_bed(wildcards, 'excluded')
	output:
		tmp = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.filter3.txt.tmp"),
		kept = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.filter3.txt"),
		excluded = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.removed3.txt"),
	resources:
		mem_mb=lambda wildcards, attempt, input: int(resources_list_with_min_and_max(input, attempt, 1.5)),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
	container:
		"docker://szsctt/bedtools:1"
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		cp {input.filt} {output.tmp}
 
		ARRAY=( {input.beds} )
		for bed in "${{ARRAY[@]}}"; do
			echo "only keepint integrations that intersect with $bed"
			head -n1 {input.filt} > {output.kept}
			bedtools intersect -u -a {output.tmp} -b $bed >> {output.kept}
			bedtools intersect -v -a {output.tmp} -b $bed >> {output.excluded}
			cp {output.kept} {output.tmp}
		done
		"""

def get_post_final(wildcards, file_type):
	assert file_type in {'kept', 'excluded'}
	exclude_beds = len(get_value_from_df(wildcards, 'bed_exclude'))
	include_beds = len(get_value_from_df(wildcards, 'bed_exclude'))
	
	# if there aren't any bed files to use for excluding
	if exclude_beds == 0 and include_beds == 0:
		if file_type == 'kept':
			return rules.post_filter.output.kept
		else:
			return rules.post_filter.output.excluded
	# if there are, then use files after excluding
	elif exclude_beds != 0:
		if file_type == 'kept':
			return rules.exclude_bed.output.kept
		else:
			return rules.exclude_bed.output.excluded
	else:
		if file_type == 'kept':
			return rules.include_bed.output.kept
		else:
			return rules.include_bed.output.excluded		


rule post_final:
	input:
		kept = lambda wildcards: get_post_final(wildcards, 'kept'),
		excluded = lambda wildcards: get_post_final(wildcards, 'excluded'),
	output:
		kept = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.post.txt",
		excluded = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.filter_fail.txt"
	resources:
		mem_mb=lambda wildcards, attempt, input: int(resources_list_with_min_and_max(input, attempt, 1.5)),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
	shell:
		"""
		mv {input.kept} {output.kept}
		mv {input.excluded} {output.excluded}
		"""

rule separate_unique_locations:
	input:
		kept = rules.post_final.output.kept
	output:
		unique = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.post.unique.txt",
		host_ambig = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.post.host_ambig.txt",
		virus_ambig = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.post.virus_ambig.txt",
		both_ambig = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.post.both_ambig.txt",
		tmp_both = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.post.both_ambig.txt.tmp")
	params:
		mapq = lambda wildcards: get_value_from_df(wildcards, 'mapq_thresh')
	container:
		"docker://szsctt/simvi:1"
	conda: "../envs/filter.yml"
	resources:
		mem_mb=lambda wildcards, attempt, input: int(resources_list_with_min_and_max(input, attempt, 1.5)),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
	shell:
		"""
		# get uniquely localised integrations as kept
		python3 scripts/filter.py \
		-i {input.kept} \
		-k {output.unique} \
		-e {output.both_ambig} \
		-c 'HostMapQ > {params.mapq} and ViralMapQ > {params.mapq} and AltLocs == None'
		
		# get integrations that are only ambiguous in host (but not in virus)
		python3 scripts/filter.py \
		-i {output.both_ambig} \
		-k {output.host_ambig} \
		-e {output.tmp_both} \
		-c 'ViralMapQ > {params.mapq} and AltLocs == uniqueVirus'	
			
		# get integrations that are only ambiguous in virus (but not in host)
		python3 scripts/filter.py \
		-i {output.tmp_both} \
		-k {output.virus_ambig} \
		-e {output.both_ambig} \
		-c 'HostMapQ > {params.mapq} and AltLocs == uniqueHost'		
		"""
		
rule merged_bed:
	input:
		txt = rules.separate_unique_locations.output.unique
	output:
		merged = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.post.unique.merged.txt"
	params:
		method = lambda wildcards: get_value_from_df(wildcards, 'merge_method'),
		n = lambda wildcards: int(get_value_from_df(wildcards, 'merge_n_min')),
		
	container:
		"docker://szsctt/bedtools:1"
	resources:
		mem_mb=lambda wildcards, attempt, input: int(resources_list_with_min_and_max(input, attempt, 1.5, 1000)),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
	threads: 1
	shell:
		"""
		python3 scripts/merge.py -i {input.txt} -o {output.merged} -c {params.method} -n {params.n}
		"""


rule summarise:
	input: 
		merged_beds = lambda wildcards: expand(strip_wildcard_constraints(rules.merged_bed.output.merged), zip,
							samp = toDo.loc[toDo['dataset'] == wildcards.dset,'sample'],
							host = toDo.loc[toDo['dataset'] == wildcards.dset,'host'],
							virus = toDo.loc[toDo['dataset'] == wildcards.dset,'virus'],
							allow_missing = True
					)
	output:
		"{outpath}/summary/{dset}.xlsx",
		"{outpath}/summary/{dset}_annotated.xlsx"
	conda:
		"../envs/rscripts.yml"
	container:
		"docker://szsctt/rscripts:4"
	params:
		outdir = lambda wildcards, output: path.dirname(output[0]),
		host = lambda wildcards: set(toDo.loc[toDo['dataset'] == wildcards.dset,'host']).pop(),
		virus = lambda wildcards: set(toDo.loc[toDo['dataset'] == wildcards.dset,'virus']).pop()
	resources:
		mem_mb=lambda wildcards, attempt, input: resources_list_with_min_and_max(input, attempt, 3, 1000),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
	threads: 1
	shell:
		"Rscript scripts/summarise_ints.R {params.host} {params.virus} {input} {params.outdir}"

rule ucsc_bed:
	input:
		merged_beds = lambda wildcards: expand(strip_wildcard_constraints(rules.merged_bed.output.merged), zip,
							samp = toDo.loc[toDo['dataset'] == wildcards.dset,'sample'],
							host = toDo.loc[toDo['dataset'] == wildcards.dset,'host'],
							virus = toDo.loc[toDo['dataset'] == wildcards.dset,'virus'],
							allow_missing = True
					)
	output:
		"{outpath}/summary/ucsc_bed/{dset}.post.bed"
	params:
		outdir = lambda wildcards, output: f"{os.path.dirname(output[0])}/{wildcards.dset}",
		host = lambda wildcards: set(toDo.loc[toDo['dataset'] == wildcards.dset,'host']).pop(),
		virus = lambda wildcards: set(toDo.loc[toDo['dataset'] == wildcards.dset,'virus']).pop()
	conda:
		"../envs/rscripts.yml"
	container:
		"docker://szsctt/rscripts:4"
	resources:
		mem_mb=lambda wildcards, attempt, input: resources_list_with_min_and_max(input, attempt, 3, 1000),
		time = lambda wildcards, attempt: (30, 120, 1440, 10080)[attempt - 1],
	threads: 1
	shell:
		"""
		Rscript scripts/writeBed.R {params.host} {params.virus} {input} {params.outdir}
		bash -e scripts/format_ucsc.sh {params.outdir}
		"""

