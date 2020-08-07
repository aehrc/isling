
rule run_soft:
	input:
		host = lambda wildcards: get_sam(wildcards, "combined", "host"),
		virus = lambda wildcards: get_sam(wildcards, "combined", "virus")
	output:
		soft = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.soft.txt"),
	container:
		"docker://ubuntu:18.04"	
	shell:
		"""
		perl -I. ./softClip.pl --viral {input.virus} --human {input.host} --output {output.soft} --tol 3
		"""
		
rule run_short:
	input:
		host = lambda wildcards: get_sam(wildcards, "combined", "host"),
		virus = lambda wildcards: get_sam(wildcards, "combined", "virus"),
	output:
		short = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.short.txt"),
	container:
		"docker://ubuntu:18.04"
	shell:
		"""
		perl -I. ./short.pl --viral {input.virus} --human {input.host} --output {output.short} --tol 3
		"""
		
rule run_discordant:
	input:
		host = lambda wildcards: get_sam(wildcards, "paired", "host"),
		virus = lambda wildcards: get_sam(wildcards, "paired", "virus"),
	output:
		discord = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.discordant.txt"),
	container:
		"docker://ubuntu:18.04"
	shell:
		"""
		perl -I. ./discordant.pl --viral {input.virus} --human {input.host} --output {output.discord} --tol 3
		"""

rule combine_ints:
	input:
		soft = rules.run_soft.output,
		short = rules.run_short.output,
		discordant = rules.run_discordant.output
	output:
		temp =  temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.txt.tmp"),
		all = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.txt"
	container:
		"docker://ubuntu:18.04"
	shell:
		"""
		awk 'FNR>1 || NR==1' {input} > {output.all}
		sort -n -k1,1 -k2,2n {output.all} > {output.temp}
		cp {output.temp} {output.all}
		"""

