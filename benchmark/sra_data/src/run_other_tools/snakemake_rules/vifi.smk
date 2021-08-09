#####################################################
####################### vifi ########################
#####################################################

rule host_virus_index:
	input:
		virus = lambda wildcards: ref_names[wildcards.virus],
		host = lambda wildcards: ref_names[wildcards.host]
	output:
		fa = "{outpath}/vifi_refs/data/{virus}/{host}_{virus}.fas",
		idx = multiext("{outpath}/vifi_refs/data/{virus}/{host}_{virus}.fas", ".amb", ".ann", ".bwt", ".pac", ".sa")
	container:
		"docker://szsctt/vifi:1"
	resources:
		mem_mb= lambda wildcards, attempt, input: resources_list_with_min_and_max((input.host, input.virus), attempt, 2, 1000),
		time = lambda wildcards, attempt: ('2:00:00', '24:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1
	shell:
		"""
		cat {input.virus} {input.host} > {output.fa}
		bwa index {output.fa}
		"""
		
rule vifi_faidx:
	input:
		fa = lambda wildcards: ref_names[wildcards.host]
	output:
		fai = "{outpath}/vifi_refs/data_repo/{host}/{host}.fa.fai"
	params:
		fai = lambda wildcards, input: input.fa + ".fai"
	container:
		"docker://szsctt/vifi:1"
	shell:
		"""
		samtools faidx {input.fa}
		mv {params.fai} {output.fai}
		"""
		
rule vifi_data_repo:
	input:
		fa = lambda wildcards: ref_names[wildcards.host],
		fai = rules.vifi_faidx.output.fai,
		mappability = lambda wildcards: get_vifi_resource(wildcards, analysis_df, 'host_mappability'),
		mappability_exclude = lambda wildcards: get_vifi_resource(wildcards, analysis_df, 'host_mappability_exclude'),
		genes = lambda wildcards: get_vifi_resource(wildcards, analysis_df, 'host_genes'),
		exons = lambda wildcards: get_vifi_resource(wildcards, analysis_df, 'host_exons'),
		oncogenes = lambda wildcards: get_vifi_resource(wildcards, analysis_df, 'host_oncogenes'),
		centromeres = lambda wildcards: get_vifi_resource(wildcards, analysis_df, 'host_centromeres'),
		conserved = lambda wildcards: get_vifi_resource(wildcards, analysis_df, 'host_conserved_regions'),
		segdup = lambda wildcards: get_vifi_resource(wildcards, analysis_df, 'host_segdup')	
	output:
		chromosomes = "{outpath}/vifi_refs/data_repo/{host}/{host}_chromosome-list.txt",
		file_list = "{outpath}/vifi_refs/data_repo/{host}/file_list.txt",
		fa = "{outpath}/vifi_refs/data_repo/{host}/{host}.fa",
		mappability = "{outpath}/vifi_refs/data_repo/{host}/{host}.mappability.bed",
		mappability_exclude = "{outpath}/vifi_refs/data_repo/{host}/{host}.mappability-exclude.bed",
		genes = "{outpath}/vifi_refs/data_repo/{host}/{host}.genes.gff",
		exons = "{outpath}/vifi_refs/data_repo/{host}/{host}.exons.gff"	,
		oncogenes = "{outpath}/vifi_refs/data_repo/{host}/{host}.oncogenes.gff",
		centromeres = "{outpath}/vifi_refs/data_repo/{host}/{host}.centromeres.bed",
		conserved = "{outpath}/vifi_refs/data_repo/{host}/{host}.conserved.bed",
		segdup = "{outpath}/vifi_refs/data_repo/{host}/{host}.segdup.bed"
	shell:
		"""
		cp {input.fa} {output.fa}
		cp {input.mappability} {output.mappability}
		cp {input.mappability_exclude} {output.mappability_exclude}
		cp {input.genes} {output.genes}
		cp {input.exons} {output.exons}
		cp {input.oncogenes} {output.oncogenes}
		cp {input.centromeres} {output.centromeres}
		cp {input.conserved} {output.conserved}
		cp {input.segdup} {output.segdup}
		
		awk 'BEGIN {{ ORS = " " }} {{a[$1]}} END {{for (i in a) print i}}' \
				{input.fai} > {output.chromosomes}
				
		echo "fa_file 		                  $(basename {output.fa})" > {output.file_list}
		echo "chrLen_file 		              $(basename {input.fai})" >> {output.file_list}
  	echo "duke35_filename 		          $(basename {output.mappability})" >> {output.file_list}
  	echo "mapability_exclude_filename   $(basename {output.mappability_exclude})" >> {output.file_list}
  	echo "gene_filename 		            $(basename {output.genes})" >> {output.file_list}
  	echo "exon_file 		                $(basename {output.exons})" >> {output.file_list}
  	echo "oncogene_filename 		        $(basename {output.oncogenes})" >> {output.file_list}
  	echo "centromere_filename 		      $(basename {output.centromeres})" >> {output.file_list}
  	echo "conserved_regions_filename 		$(basename {output.conserved})" >> {output.file_list}
  	echo "segdup_filename 		          $(basename {output.segdup})" >> {output.file_list}		
		"""
	
def vifi_hosts(wildcards):
	"""get all the names of the hosts used for vifi"""
	subset = analysis_df[(analysis_df['tool'] == 'vifi') & (analysis_df['outdir'] == wildcards.outpath)]
	return set(subset['host'])
	
rule data_repo_host_list:
	input: 
		lambda wildcards: [f"{wildcards.outpath}/vifi_refs/data_repo/{host}/file_list.txt" for host in vifi_hosts(wildcards)]
	output: "{outpath}/vifi_refs/data_repo/reference.txt"
	run:
		shell("rm -f {output}")
		shell("touch {output}")	
		hosts = vifi_hosts(wildcards)	
		for host in hosts:
			shell("echo {host} >> {output}")
		
rule vifi:
	input:
		rules.vifi_data_repo.output.chromosomes,
		rules.vifi_data_repo.output.file_list,				
		rules.vifi_data_repo.output.mappability,	
		rules.vifi_data_repo.output.mappability_exclude,
		rules.vifi_data_repo.output.genes,		
		rules.vifi_data_repo.output.exons,		
		rules.vifi_data_repo.output.oncogenes,
		rules.vifi_data_repo.output.centromeres,		
		rules.vifi_data_repo.output.conserved,		
		rules.vifi_data_repo.output.segdup,		
		rules.vifi_faidx.output.fai,
		fa = rules.vifi_data_repo.output.fa,
		idx = rules.host_virus_index.output.idx,
		host_list = rules.data_repo_host_list.output[0],
		fastq1 = lambda wildcards: get_reads(wildcards, analysis_df, rules, 1),
		fastq2 = lambda wildcards: get_reads(wildcards, analysis_df, rules, 2),
		chrom_list = rules.vifi_data_repo.output.chromosomes,
	output:
		clusters = "{outpath}/{dset}/vifi/{samp}.{host}.{virus}/output.clusters.txt",
		range = "{outpath}/{dset}/vifi/{samp}.{host}.{virus}/output.clusters.txt.range"
	params:
		reference_repo = lambda wildcards, input: os.path.dirname(input.fa),
		aa_data_repo = lambda wildcards, input: os.path.dirname(input.host_list),
		outdir = lambda wildcards, output: os.path.dirname(output.clusters),
		reference = lambda wildcards, input: os.path.splitext(input.idx[0])[0]
	wildcard_constraints:
		analysis_condition = "vifi\d+"
	container:
		"docker://szsctt/vifi:1"
	resources:
		mem_mb= lambda wildcards, attempt, input: resources_list_with_min_and_max((input.idx), attempt),
		time = lambda wildcards, attempt: ('2:00:00', '24:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1
	threads: 8
	shell:
		"""
		export AA_DATA_REPO=$(realpath {params.aa_data_repo})
		export REFERENCE_REPO=$(realpath {params.reference_repo})
		
		python $VIFI_DIR/scripts/run_vifi.py \
		-c {threads} \
		-f $(realpath {input.fastq1}) -r $(realpath {input.fastq2}) \
		--reference $(realpath {params.reference}) \
		-v {wildcards.virus} \
		-o $(realpath {params.outdir}) \
		-d True \
		-C $(realpath {input.chrom_list})
		"""			
