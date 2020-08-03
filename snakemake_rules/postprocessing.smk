
rule sortbed:
	input:
		TOSORT
	output:
		SORTED
	run:
		for i in range(len(TOSORT)):
			print(f"sorting file {TOSORT[i]} into file {SORTED[i]}, file extension {path.splitext(TOSORT[i])[1]}")
			
			# if file is a bed file
			if (path.splitext(TOSORT[i])[1] == ".bed"):
				body = f"sort -k1,1 -k2,2n {TOSORT[i]} > {SORTED[i]}"
				shell(body)
				print(body)
			# if file is a gtf file
			elif (path.splitext(TOSORT[i])[1] == ".gtf"):
				
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
		"{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.txt",
		rules.sortbed.output
	output:
		"{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.post.txt"
	conda:
		"../envs.rscripts.yml"
	container:
		"docker://szsctt/rscripts:2"
	params:
		lambda wildcards, input: f"{input[0]} {POSTARGS[wildcards.dset]}"
	shell:
		"""
		Rscript post/postprocess.R {params}
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
	conda:
		"../envs.rscripts.yml"
	container:
		"docker://szsctt/rscripts:2"
	params:
		outdir = lambda wildcards, output: path.dirname(output[0])
	shell:
		"Rscript summarise_ints.R {input} {params.outdir}"

rule write_bed:
	input:
		

rule ucsc_bed:
	input:
		lambda wildcards: [f"{wildcards.outpath}/{wildcards.dset}/ints/{samp}.{host}.{virus}.integrations.post.txt" for samp, host, virus
			in zip(toDo.loc[toDo['dataset'] == wildcards.dset,'sample'], 
					toDo.loc[toDo['dataset'] == wildcards.dset,'host'], 
					toDo.loc[toDo['dataset'] == wildcards.dset,'virus'])]
	output:
		"{outpath}/summary/ucsc_bed/{dset}.post.bed"
	params:
		outdir = lambda wildcards, output: f"{path.dirname(output[0])}/{wildcards.dset}"
	conda:
		"../envs.rscripts.yml"
	container:
		"docker://szsctt/rscripts:2"
	shell:
		"""
		Rscript writeBed.R {input} {params.outdir}
		bash -e format_ucsc.sh {params.outdir}
		mv {params.outdir}/*bed {params.outdir}/..
		rmdir {params.outdir}
		"""
