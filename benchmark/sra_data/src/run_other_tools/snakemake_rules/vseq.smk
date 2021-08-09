#####################################################
################### VSeq-Toolkit ####################
#####################################################

rule bwa_index:
	input:
		genome = lambda wildcards: ref_names[wildcards.genome],
	output:
		idx = multiext("{outpath}/refs/bwa/{genome}/{genome}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa")
	params:
		prefix = lambda wildcards, output: os.path.splitext(output[0])[0]
	container:
		"docker://szsctt/vseq:1"
	resources:
		mem_mb= lambda wildcards, attempt, input: resources_list_with_min_and_max((input.genome, ), attempt, 2, 1000),
		time = lambda wildcards, attempt: ('2:00:00', '24:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1	
	shell:
		"""
		 $VSeqToolkit/thirdPartyTools/bwa index -p {params.prefix} {input.genome}
		"""

rule link_virus:
	input:
		genome = lambda wildcards: ref_names[wildcards.genome],
	output:
		genome = "{outpath}/refs/bwa/{genome}/{genome}.fa"
	container:
		"docker://szsctt/vseq:1"
	shell:
		"""
		cp $(realpath {input.genome}) $(realpath {output.genome})
		"""
	
def format_sed_filename(to_replace, start, val):
	return f"-e 's#{to_replace}#{start}{val}#g'"

def format_sed_command(wildcards, analysis_df, to_replace, start, column):
	val =  analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', column)
	return format_sed_filename(to_replace, start, val)
	
def format_sed_command_int(wildcards, analysis_df, to_replace, start, column):
	val =  int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', column))
	return format_sed_filename(to_replace, start, val)
		
rule vseq_toolkit_config_template:
	input:
		combined_idx = rules.host_virus_index_seeksv.output.idx,
		vec_idx = multiext("{outpath}/refs/bwa/{virus}/{virus}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
		fastq1 = lambda wildcards: get_input_reads(wildcards, analysis_df, 1),
		fastq2 = lambda wildcards: get_input_reads(wildcards, analysis_df, 2),
	output:
		config = "{outpath}/{dset}/vseq_toolkit/{samp}.{host}.{virus}.config.txt"
	params:
		fastq1 = lambda wildcards, input: format_sed_filename(
								"file1= $VSeqToolkit/testDir/testData/testDataCombined.R1.fastq.gz", 
								"file1= ", os.path.abspath(input.fastq1)
								),
		fastq2 = lambda wildcards, input: format_sed_filename(
								"file2= $VSeqToolkit/testDir/testData/testDataCombined.R2.fastq.gz", 
								"file2= ", os.path.abspath(input.fastq2)
								),	
		out_dir = lambda wildcards, output: format_sed_filename(
								"outDir= $VSeqToolkit/testDir/testResultsCheck/", 
								"outDir= ", f"{os.path.abspath(os.path.splitext(os.path.splitext(output.config)[0])[0])}/",
								),
		vec_idx = lambda wildcards, input: format_sed_filename(
								"vecRef= $VSeqToolkit/testDir/testReferenceIndex/vector1.fa", 
								"vecRef= ", os.path.abspath(os.path.splitext(input.vec_idx[0])[0]),
								),
		combined_idx = lambda wildcards, input: format_sed_filename(
								"vecGenRef= $VSeqToolkit/testDir/testReferenceIndex/hg38chr22Vector1.fa", 
								"vecGenRef= ", os.path.abspath(os.path.splitext(input.combined_idx[0])[0]),
								),
		annoTable = lambda wildcards: format_sed_command(wildcards, analysis_df, 
								"annoTable= $VSeqToolkit/testDir/testReferenceIndex/refSeqUCSCTablehg38.txt", 
								"annoTable= ", 'host_table'),

		
		adapter1 = lambda wildcards: format_sed_command(wildcards, analysis_df, 
																										"adapter1=GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", 
																											"adapter1=", 'adapter_1'),
		adapter2 = lambda wildcards: format_sed_command(wildcards, analysis_df, 
																		"adapter2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT", 
																		"adapter2=", 'adapter_2'),
		qual = lambda wildcards: format_sed_command_int(wildcards, analysis_df, 
																								"qua=20","qua=", 'qual'),
		lenPer = lambda wildcards: format_sed_command_int(wildcards, analysis_df, 
																									"lenPer=50","lenPer=", 'lenPer'),
		mode = lambda wildcards: format_sed_command(wildcards, analysis_df, 
																									"mode=default","mode=", 'mode'),
		contAna = lambda wildcards: format_sed_filename("contAna=true", "contAna=", "false"),
		vecVecFusion = lambda wildcards: format_sed_command(wildcards, analysis_df, 
																											"vecVecFusion=true","vecVecFusion=",
																												 'vecVecFusion'),
		stringencyVec = lambda wildcards: format_sed_command(wildcards, analysis_df, 
																													"stringencyVec=low","stringencyVec=",
																													 'stringencyVec'),
		UMthresholdVec = lambda wildcards: format_sed_command(wildcards, analysis_df, 
																														"UMthresholdVec=0.95","UMthresholdVec=",
																														 'UMthresholdVec'),		
		minMapSpanVec = lambda wildcards: format_sed_command_int(wildcards, analysis_df, 
																														"minMapSpanVec=20","minMapSpanVec=",
																														'minMapSpanVec'),	
		distVecVec = lambda wildcards: format_sed_command_int(wildcards, analysis_df, 
																												"distVecVec=10","distVecVec=",
																												'distVecVec'),	
		opVecVec = lambda wildcards: format_sed_command_int(wildcards, analysis_df, 
																											"opVecVec=5","opVecVec=",
																											'opVecVec'),
		idenVecVec = lambda wildcards: format_sed_command_int(wildcards, analysis_df, 
																											"idenVecVec=95","idenVecVec=",
																											'idenVecVec'),		
		stringencyVecGen = lambda wildcards: format_sed_command(wildcards, analysis_df, 
																														"stringencyVecGen=low","stringencyVecGen=",
																														 'stringencyVecGen'),			
		UMthresholdVecGen = lambda wildcards: format_sed_command(wildcards, analysis_df, 
																															"UMthresholdVecGen=0.95",
																															"UMthresholdVecGen=",
																														 	'UMthresholdVecGen'),	
		minMapSpanVecGen = lambda wildcards: format_sed_command_int(wildcards, analysis_df, 
																														"minMapSpanVecGen=20",
																														"minMapSpanVecGen=",
																														'minMapSpanVecGen'),	
		distVecGen = lambda wildcards: format_sed_command_int(wildcards, analysis_df, 
																											"distVecGen=10","distVecGen=",
																											'distVecGen'),	
		opVecGen = lambda wildcards: format_sed_command_int(wildcards, analysis_df, 
																										"opVecGen=5","opVecGen=",
																										'opVecGen'),			
		idenVecGen = lambda wildcards: format_sed_command_int(wildcards, analysis_df, 
																												"idenVecGen=95","idenVecGen=",
																												'idenVecGen'),			
		clusterRange = lambda wildcards: format_sed_command_int(wildcards, analysis_df, 
																													"clusterRange=3","clusterRange=",
																													'clusterRange'),	 
	container:
		"docker://szsctt/vseq:1"
	shell:
		"""
		cp /var/work/VSeq-Toolkit/config.test.txt {output.config}
		sed -i {params} {output.config}
		"""
		
rule vseq_toolkit:
	input:
		combined_idx = rules.host_virus_index_seeksv.output.idx,
		combined_genome = rules.host_virus_index_seeksv.output.fa,
		vec_idx = multiext("{outpath}/refs/bwa/{virus}/{virus}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
		vec_genome = "{outpath}/refs/bwa/{virus}/{virus}.fa",
		fastq1 = lambda wildcards: get_input_reads(wildcards, analysis_df, 1),
		fastq2 = lambda wildcards: get_input_reads(wildcards, analysis_df, 2),
		config = rules.vseq_toolkit_config_template.output.config
	output:
		clust = "{outpath}/{dset}/vseq_toolkit/{samp}.{host}.{virus}/ISGenomeVector.csv",
		nonUnique = "{outpath}/{dset}/vseq_toolkit/{samp}.{host}.{virus}/ISGenomeVector.NonUniqueGenome.csv",
		unique = "{outpath}/{dset}/vseq_toolkit/{samp}.{host}.{virus}/ISGenomeVector.UniqueGenome.csv",
		unclust = "{outpath}/{dset}/vseq_toolkit/{samp}.{host}.{virus}/ISGenomeVector.Unclustered.csv"
	container:
		"docker://szsctt/vseq:1"
	threads: 10 # hard-coded into scripts that run bwa-mem?
	resources:
		mem_mb= lambda wildcards, attempt, input: resources_list_with_min_and_max((input.combined_genome, input.vec_genome), attempt, 20000),
		time = lambda wildcards, attempt: ('30:00', '2:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1	
	shell:
		"""
		perl -I $VSeqToolkit/scripts/ $VSeqToolkit/scripts/VSeq-TK.pl -c $(realpath {input.config})
		
		if [ ! -e {output.clust} ]; then
			touch {output.clust}
			touch {output.nonUnique}
			touch {output.unique}
			touch {output.unclust}
		fi
		"""		
		
				
