import os

#####################################################
###################### verse ########################
#####################################################

# verse uses a different version of bowtie2 than polyidus
# in theory they could share the same indexes
# but just to be safe, index again
rule bwt2_verse:
	input:
		fasta = lambda wildcards: ref_names[wildcards.genome]
	output:
		multiext("{outpath}/verse_references/{genome}/{genome}", 
							".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"
						)
	container:
		"docker://szsctt/verse:1"
	params:
		prefix = lambda wildcards, output: os.path.splitext(os.path.splitext(output[0])[0])[0]
	resources:
		mem_mb= lambda wildcards, attempt, input: resources_list_with_min_and_max((input.fasta, ), attempt, 2, 1000),
		time = lambda wildcards, attempt: ('2:00:00', '24:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1
	shell:
		"bowtie2-build {input} {params.prefix}"
		
rule blastplus_verse:
	input:
		fasta = lambda wildcards: ref_names[wildcards.genome]
	output:
		multiext("{outpath}/verse_references/{genome}/{genome}", 
							".nhr", ".nin", ".nsq"
						)
	container:
		"docker://szsctt/verse:1"
	params:
		prefix = lambda wildcards, output: os.path.splitext(output[0])[0]
	resources:
		mem_mb= lambda wildcards, attempt, input: resources_list_with_min_and_max((input.fasta, ), attempt),
		time = lambda wildcards, attempt: ('2:00:00', '24:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1
	shell:
		"makeblastdb -in {input} -dbtype nucl -out {params.prefix}"
		
rule copy_reference_verse:
	input:
		fasta = lambda wildcards: ref_names[wildcards.genome],
	output:
		fasta = "{outpath}/verse_references/{genome}/{genome}.fa",
	shell:
		"""
		cp $(realpath {input.fasta}) $(realpath {output.fasta})
		"""

rule verse:
	input:
		fastq1 = lambda wildcards: get_reads(wildcards, analysis_df, rules, 1),
		fastq2 = lambda wildcards: get_reads(wildcards, analysis_df, rules, 2),
		virus_fasta = "{outpath}/verse_references/{virus}/{virus}.fa",
		host_fasta = "{outpath}/verse_references/{host}/{host}.fa",
		host_idx = multiext("{outpath}/verse_references/{host}/{host}", 
												".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"
												),
		virus_blast = multiext("{outpath}/verse_references/{virus}/{virus}", 
														".nhr", ".nin", ".nsq"),
		host_blast = multiext("{outpath}/verse_references/{host}/{host}", 
														".nhr", ".nin", ".nsq"),
	output:
		output = "{outpath}/{dset}/verse/{samp}.{host}.{virus}/integration-sites.txt",
		config = "{outpath}/{dset}/verse/{samp}.{host}.{virus}/config.txt",	
	params:
		workdir = lambda wildcards, output: os.path.dirname(os.path.realpath(output.output)),
		currdir = lambda wildcards: os.getcwd(),
		config = lambda wildcards, output: os.path.basename(output.config),
		host_bowtie = lambda wildcards, input: os.path.splitext(os.path.splitext(os.path.realpath(input.host_idx[0]))[0])[0],
		virus_blast = lambda wildcards, input: os.path.splitext(os.path.realpath(input.virus_blast[0]))[0],
		host_blast = lambda wildcards, input: os.path.splitext(os.path.realpath(input.host_blast[0]))[0],
		virus_fa = lambda wilcards, input: os.path.realpath(input.virus_fasta),
		detection_mode = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'verse', 'detection_mode'),
		flank_region_size = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'verse', 'flank_region_size')),
		sensitivity_level = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'verse', 'sensitivity_level')),
		min_contig_length = lambda wildcards: int(analysis_df_tool_value(wildcards,analysis_df, 'verse', 'min_contig_length')),
		blastn_evalue_thrd = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'verse', 'blastn_evalue_thrd'),
		similarity_thrd = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'verse', 'similarity_thrd'),
		chop_read_length = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'verse', 'chop_read_length')),
		minIdentity = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'verse', 'minIdentity')),	
	container:
		"docker://szsctt/verse:1"
	resources:
		mem_mb= 3000000,
		time= '7-00:00:00',
		nodes = 1
	log:
		"{outpath}/logs/{dset}_verse_{samp}_{host}_{virus}.log"
	threads: verse_threads
	shell:
		"""
		mkdir -p {params.workdir}
		touch {output.config}
		
		# write config file
		echo "fastq1 = $(readlink -f {input.fastq1})" > {output.config}
		echo "fastq2 = $(readlink -f {input.fastq2})" >> {output.config}
		echo 'detect_integration = yes' >> {output.config}
		echo 'detect_mutation = no' >> {output.config}	
		echo 'thread_no = {verse_threads}' >> {output.config}
  		echo 'blastn_bin = /usr/bin/blastn' >> {output.config}
		echo 'bowtie_bin = /usr/bin/bowtie2' >> {output.config}
		echo 'bwa_bin = /usr/bin/bwa' >> {output.config}
		echo 'trinity_script = /trinityrnaseq_r2013-02-16/Trinity.pl' >> {output.config}
		echo 'SVDetect_dir = /SVDetect_r0.8' >> {output.config}
		echo 'virus_database = {params.virus_fa}' >> {output.config}
		echo 'bowtie_index_human = {params.host_bowtie}' >> {output.config}
		echo 'blastn_index_human = {params.host_blast}' >> {output.config}
		echo 'blastn_index_virus = {params.virus_blast}' >> {output.config}
		echo 'detection_mode     = {params.detection_mode}' >> {output.config}
		echo 'flank_region_size  = {params.flank_region_size}' >> {output.config}
		echo 'sensitivity_level  = {params.sensitivity_level}' >> {output.config}
		echo 'min_contig_length = {params.min_contig_length}' >> {output.config}
		echo 'blastn_evalue_thrd = {params.blastn_evalue_thrd}' >> {output.config}
		echo 'similarity_thrd = {params.similarity_thrd}' >> {output.config}
		echo 'chop_read_length = {params.chop_read_length}' >> {output.config}
		echo 'minIdentity = {params.minIdentity}' >> {output.config}

		cd {params.workdir}
			
		perl /var/work/VirusFinder2.0/VirusFinder.pl -c $(basename {output.config})
		
		cd {params.currdir}
		# virusFinder only makes an output file if it finds integration sites
		# but snakemake always requires the same output files
		if [ ! -e {output.output} ]; then
			touch {output.output}
		fi
		"""
		

