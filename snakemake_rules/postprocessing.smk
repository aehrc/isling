rule filter:
	input:
		ints = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.txt",
	output:
		kept = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.p.txt"),
		excluded = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.re.txt"),
	params:
		filterstring = lambda wildcards: get_value_from_df(wildcards, 'filter')
	container:
		"docker://szsctt/simvi:1"
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
		filt = rules.filter.output.kept,
		excluded = rules.filter.output.excluded
	output:
		tmp = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.po.txt.tmp"),
		kept = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.po.txt",
		excluded = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.remo.txt",
	container:
		"docker://szsctt/bedtools:1"
	conda:
		"envs/bedtools.yml"
	shell:
		"""
		cp {input.filt} {output.tmp}
		cp {input.excluded} {output.excluded}
		
		if [ -z "{input.beds}" ]; then
			cp {input.filt} {output.kept}
			cp {input.excluded} {output.excluded}
			touch {output.tmp}
		else
			ARRAY=( {input.beds} )
			for bed in "${{ARRAY[@]}}"; do
				echo "excluding integrations that intersect with $bed"
				head -n1 {input.filt} > {output.kept}
				bedtools intersect -v -a {output.tmp} -b $bed >> {output.kept}
				bedtools intersect -u -a {output.tmp} -b $bed >> {output.excluded}
				cp {output.kept} {output.tmp}
			done
		fi
		"""

rule include_bed:
	input:
		beds = lambda wildcards: [os.path.splitext(i)[0] + '.sorted.bed' for i in get_value_from_df(wildcards, 'bed_include')],	
		filt = rules.exclude_bed.output.kept,
		excluded = rules.exclude_bed.output.excluded
	output:
		tmp = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.post.txt.tmp"),
		kept = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.post.txt",
		excluded = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.removed.txt",
	container:
		"docker://szsctt/bedtools:1"
	conda:
		"envs/bedtools.yml"
	shell:
		"""
		cp {input.filt} {output.tmp}
		cp {input.excluded} {output.excluded}
		
		if [ -z "{input.beds}" ]; then
			cp {input.filt} {output.kept}
			cp {input.excluded} {output.excluded}
			touch {output.tmp}
		else
			ARRAY=( {input.beds} )
			for bed in "${{ARRAY[@]}}"; do
				echo "only keepint integrations that intersect with $bed"
				head -n1 {input.filt} > {output.kept}
				bedtools intersect -u -a {output.tmp} -b $bed >> {output.kept}
				bedtools intersect -v -a {output.tmp} -b $bed >> {output.excluded}
				cp {output.kept} {output.tmp}
			done
		fi
		"""

rule summarise:
	input:
		lambda wildcards: [f"{{outpath}}/{{dset}}/ints/{samp}.{host}.{virus}.integrations.post.txt" for samp, host, virus
			in zip(toDo.loc[toDo['dataset'] == wildcards.dset,'sample'], 
					toDo.loc[toDo['dataset'] == wildcards.dset,'host'], 
					toDo.loc[toDo['dataset'] == wildcards.dset,'virus'])]
	output:
		"{outpath}/summary/{dset}.xlsx",
		"{outpath}/summary/{dset}_annotated.xlsx"
#	group: "post"
	conda:
		"../envs/rscripts.yml"
	container:
		"docker://szsctt/rscripts:4"
	params:
		outdir = lambda wildcards, output: path.dirname(output[0])
	resources:
		mem_mb=lambda wildcards, attempt, input: resources_list_with_min_and_max(input, attempt, 3, 1000)
	threads: 1
	shell:
		"Rscript scripts/summarise_ints.R {input} {params.outdir}"


rule ucsc_bed:
	input:
		lambda wildcards: [f"{wildcards.outpath}/{wildcards.dset}/ints/{samp}.{host}.{virus}.integrations.post.txt" for samp, host, virus
			in zip(toDo.loc[toDo['dataset'] == wildcards.dset,'sample'], 
					toDo.loc[toDo['dataset'] == wildcards.dset,'host'], 
					toDo.loc[toDo['dataset'] == wildcards.dset,'virus'])]
	output:
		"{outpath}/summary/ucsc_bed/{dset}.post.bed"
#	group: "post"
	params:
		outdir = lambda wildcards, output: f"{os.path.dirname(output[0])}/{wildcards.dset}"
	conda:
		"../envs/rscripts.yml"
	container:
		"docker://szsctt/rscripts:4"
	resources:
		mem_mb=lambda wildcards, attempt, input: resources_list_with_min_and_max(input, attempt, 3, 1000)
	threads: 1
	shell:
		"""
		Rscript scripts/writeBed.R {input} {params.outdir}
		bash -e scripts/format_ucsc.sh {params.outdir}
		mv {params.outdir}/*bed {params.outdir}/..
		rmdir {params.outdir}
		"""
		
rule merged_bed:
	input:
		txt = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations{post}.txt"
	output:
		sorted = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations{post}.sorted.txt"),
		merged = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations{post}.merged.txt"
#	group: "post"
	params:
		d = lambda wildcards: int(get_value_from_df(wildcards, 'merge_dist')),
		n = lambda wildcards: int(get_value_from_df(wildcards, 'merge_n_min')),
	container:
		"docker://szsctt/bedtools:1"
	resources:
		mem_mb=lambda wildcards, attempt, input: resources_list_with_min_and_max(input, attempt, 3, 1000)
	threads: 1
	shell:
		"""
		pwd
		head -n1 {input.txt} > {output.sorted}
		tail -n+2  {input.txt} | sort -k1,1 -k2,2n -k4,4 -k5,5n >> {output.sorted}
		python3 scripts/merge.py -i {output.sorted} -o {output.merged} -d {params.d} -n {params.n}
		"""

