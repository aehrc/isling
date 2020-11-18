
rule sortbed:
	input:
		TOSORT
	output:
		temp(SORTED)
	run:
		for i in range(len(TOSORT)):
			print(f"sorting file {TOSORT[i]} into file {SORTED[i]}, file extension {os.path.splitext(TOSORT[i])[1]}")
			
			# if file is a bed file
			if (os.path.splitext(TOSORT[i])[1] == ".bed"):
				body = f"sort -k1,1 -k2,2n {TOSORT[i]} > {SORTED[i]}"
				shell(body)
				print(body)
			# if file is a gtf file
			elif (os.path.splitext(TOSORT[i])[1] == ".gtf"):
				
				body = f"awk '{{{{ if ($0 !~ /^#/) {{{{ print $0 }}}} }}}}' {TOSORT[i]} | sort -k1,1 -k4,4n > {SORTED[i]}"
				shell(body)
				print(body)
			else:
				raise ValueError("only gtf files and bed files are supported")


#this rule (post) sometimes causes issues with conda. The error is usually something to do with ldpaths:

#Activating conda environment: /scratch1/sco305/intvi_cmri/intvi_pipeline/.snakemake/conda/586e76e5
#/scratch1/sco305/intvi_cmri/intvi_pipeline/.snakemake/conda/586e76e5/lib/R/bin/R: line 238: /scratch1/sco305/intvi_cmri/#intvi_pipeline/.snakemake/conda/586e76e5/lib/R/etc/ldoutpaths: No such file or directory

# however, sometimes this rule runs just fine.

# this issue is descrived here: https://github.com/conda-forge/r-base-feedstock/issues/67
# however, it doesn't appear to have been resolved.  temporarily get around this by re-attempting jobs
# that failed using command line option --restart-times, but will need to come up with a better solution for this

# I would recommend always running with singularity (--use-singularity)

rule post:
	input:
		ints = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.txt",
		sorted_beds = rules.sortbed.output
	output:
		ints = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.post.txt"
	group: "post"
	conda:
		"../envs/rscripts.yml"
	container:
		"docker://szsctt/rscripts:4"
	params:
		lambda wildcards: get_value_from_df(wildcards, 'postargs')
	shell:
		"""
		Rscript scripts/post/postprocess.R {input.ints} {params}
		"""
	
rule summarise:
	input:
		lambda wildcards: [f"{wildcards.outpath}/{wildcards.dset}/ints/{samp}.{host}.{virus}.integrations.post.txt" for samp, host, virus
			in zip(toDo.loc[toDo['dataset'] == wildcards.dset,'sample'], 
					toDo.loc[toDo['dataset'] == wildcards.dset,'host'], 
					toDo.loc[toDo['dataset'] == wildcards.dset,'virus'])]
	output:
		"{outpath}/summary/{dset}.xlsx",
		"{outpath}/summary/{dset}_annotated.xlsx"
	group: "post"
	conda:
		"../envs/rscripts.yml"
	container:
		"docker://szsctt/rscripts:4"
	params:
		outdir = lambda wildcards, output: os.path.dirname(output[0])
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
	group: "post"
	params:
		outdir = lambda wildcards, output: f"{os.path.dirname(output[0])}/{wildcards.dset}"
	conda:
		"../envs/rscripts.yml"
	container:
		"docker://szsctt/rscripts:4"
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
		bed = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations{post}.bed"),
		merged_bed = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations{post}.merged.bed"
	group: "post"
	params:
		d = lambda wildcards: f"-d {int(get_value_from_df(wildcards, 'merge_dist'))}",
		n = lambda wildcards: int(get_value_from_df(wildcards, 'merge_n_min')),
	container:
		"docker://szsctt/bedtools:1"
	shell:
		"""
		awk -F"\t" -v OFS="\t" 'BEGIN {{getline}}{{ ($9 ~ 'hv') ? dir = "+" : dir = "-"; print $1,$2,$3,dir,$21 }}' {input.txt} | sort -k1,1 -k2,2n > {output.bed}
		bedtools merge -i {output.bed} {params.d} -c 5,5 -o count,collapse | awk -F"\t" -v OFS="\t" '$4 >= {params.n}'  > {output.merged_bed}
		"""

