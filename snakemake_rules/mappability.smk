rule genmap_idx:
	input:
		fa = lambda wildcards: ref_names[wildcards.genome]
	output:
		multiext("{outpath}/references/mappability/{genome}/index.",
			 "ids.concat", "ids.limits", "info.concat", 
			 "info.limits", "lf.drp", "lf.drp.sbl", "lf.drs",
			 "lf.drv", "lf.drv.sbl", "lf.pst", "rev.lf.drp",
			 "rev.lf.drp.sbl", "ref.lf.drs", "ref.lf.drv",
			 "rev.lf.drv.sbl", "rev.lf.pst", "sa.ind", "sa.len",
			 "sa.val", "txt.concat", "txt.limits")
	params:
		dir = lambda wildcards, output: os.path.dirname(output[0])
	conda:
		"../envs/genmap.yml"
	container:
		"docker://szsctt/isling:master"
	resources:
		mem_mb = lambda wildcards, attempt, input: resources_list_with_min_and_max((input.fa,), attempt, 10, 2000),
		time =  lambda wildcards, attempt: (5, 120, 1440, 10080)[attempt - 1],
		nodes = 1
	shell:
		"genmap index -F {input.fa} -I {params.dir}"
