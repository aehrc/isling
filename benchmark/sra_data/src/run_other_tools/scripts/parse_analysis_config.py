import pdb
import itertools
import pandas as pd
from os import path
from glob import glob

# defaults
polyidus_default_aligner = ['bowtie2']
polyidus_default_trim = [False]

vifi_default_trim = [False]

seeksv_default_trim = [False]
seeksv_default_dedup = [False]

verse_default_trim = [False]
verse_default_detection_mode = ['sensitive']
verse_default_flank_region_size = [4000]
verse_default_sensitivity_level = [1]
verse_default_min_contig_length = [300]
verse_default_blastn_evalue_thrd = [0.05]
verse_default_similarity_thrd = [0.8]
verse_default_chop_read_length = [25]
verse_default_minIdentity = [80]

vseq_default_qua = [20]
vseq_default_lenPer = [50]
vseq_default_mode = ['default']
vseq_default_vecVecFusion = ['false']
vseq_default_stringencyVec = ['high']
vseq_default_UMthresholdVec = [0.95]
vseq_default_minMapSpanVec = [20]
vseq_default_distVecVec = [10]
vseq_default_opVecVec = [5]
vseq_default_idenVecVec = [0.95]
vseq_default_stringencyVecGen = ['high']
vseq_default_UMthresholdVecGen = [0.95]
vseq_default_minMapSpanVecGen = [20]
vseq_default_distVecGen = [10]
vseq_default_opVecGen = [5]
vseq_default_idenVecGen = [95]
vseq_default_clusterRange = [3]

def parse_analysis_config(config):
	
	# global parameters (defaults for this dataset) may be specified using a 'global' key
	# at the top of the config file.  These will be applied to any missing keys in the remaining 
	# datasets in the config file.  If a key is specified in both 'global' and a dataset,
	# the value specified in the dataset will be used
	
	if 'global' in config:		
		# get default (global) options
		default = config.pop('global')
		for dataset in config:
			for key in default:
				if key not in config[dataset]:
					config[dataset][key] = default[key]
					
	# get unique analysis conditions - these are combinations of the analysis parameters that
	# can be set in our pipeline (merging, de-duplicaiton, bwa mem prarameters, etc), or
	# the tool to be used in analysis

	column_names = ('experiment', 'exp', 'analysis_condition', 
									'tool', 'host', 'host_fasta',
									'virus', 'virus_fasta', 'read_folder', 
									'R1_suffix', 'R2_suffix', 'outdir', 'bam_suffix',
									'adapter_1', 'adapter_2', 'merge', 'trim', 'dedup',
									'host_mappability', 'host_mappability_exclude', 
									'host_genes', 'host_exons', 'host_oncogenes', 
									'host_centromeres', 'host_conserved_regions',
									'host_segdup', 'detection_mode', 
									'flank_region_size', 'sensitivity_level', 
									'min_contig_length', 'blastn_evalue_thrd', 
									'similarity_thrd', 
									'chop_read_length', 'minIdentity', 'qual', 'lenPer',
									'mode', 'vecVecFusion', 'stringencyVec', 'UMthresholdVec',
									'minMapSpanVec', 'distVecVec', 'opVecVec', 'idenVecVec',
									'stringencyVecGen', 'UMthresholdVecGen', 'minMapSpanVecGen', 'distVecGen', 'opVecGen',
									'idenVecGen', 'clusterRange', 'host_table'
									)		



	analysis_conditions = []
	
	for dataset in config.keys():
		
		if 'polyidus_params' in config[dataset]:
			analysis_conditions += make_polyidus_rows(config, dataset)
			
		if 'vifi_params' in config[dataset]:
			analysis_conditions += make_vifi_rows(config, dataset)
			
		if 'verse_params' in config[dataset]:
			analysis_conditions += make_verse_rows(config, dataset)

		if 'seeksv_params' in config[dataset]:
			analysis_conditions += make_seeksv_rows(config, dataset)			

		if 'vseq_toolkit_params' in config[dataset]:
			analysis_conditions += make_vseq_rows(config, dataset)		
	
	# make data frame 
	return pd.DataFrame(analysis_conditions, columns = column_names)
			
def make_verse_rows(config, dataset):
	#### parameters for verse ####
	rows = []
	
	# paramters should be lists, so that we can do combinations of all
	trim_list = get_list_with_default(config[dataset]['verse_params'], 
																		'trim', verse_default_trim, 'verse_params')	
	detection_mode_list = get_list_with_default(config[dataset]['verse_params'], 
																							'detection_mode', verse_default_detection_mode, 'verse_params')	
	flank_region_size_list = get_list_with_default(config[dataset]['verse_params'], 
																							'flank_region_size', 
																							verse_default_flank_region_size , 'verse_params')	
	sensitivity_level_list = get_list_with_default(config[dataset]['verse_params'], 
																							'sensitivity_level', 
																							verse_default_sensitivity_level, 'verse_params')		
	min_contig_length_list = get_list_with_default(config[dataset]['verse_params'], 
																							'min_contig_length', 
																							verse_default_min_contig_length, 'verse_params')
	blastn_evalue_thrd_list = get_list_with_default(config[dataset]['verse_params'], 
																							'blastn_evalue_thrd', 
																							verse_default_blastn_evalue_thrd, 'verse_params')
	similarity_thrd_list = get_list_with_default(config[dataset]['verse_params'], 
																							'similarity_thrd', 
																							verse_default_similarity_thrd, 'verse_params')
	chop_read_length_list = get_list_with_default(config[dataset]['verse_params'], 
																							'chop_read_length', 
																							verse_default_chop_read_length, 'verse_params')
	minIdentity_list = get_list_with_default(config[dataset]['verse_params'], 
																							'minIdentity', 
																							verse_default_minIdentity, 'verse_params')

 
	# each combination of these are a unique 'analysis condition' for our pipeline
	i = 0 
	for host, virus, trim, detection_mode, flank_region_size, sensitivity_level, min_contig_length, blastn_evalue_thrd, similarity_thrd, chop_read_length, minIdentity in itertools.product(
																		config[dataset]['analysis_hosts'].keys(),
																		config[dataset]['analysis_viruses'].keys(),
																		trim_list,
																		detection_mode_list,
																		flank_region_size_list,
																		sensitivity_level_list,
																		min_contig_length_list,
																		blastn_evalue_thrd_list,
																		similarity_thrd_list,
																		chop_read_length_list,
																		minIdentity_list
																		):
		
		condition = f"verse{i}"
		rows.append({
				'experiment' : dataset,
				'host' 			: host,
				'host_fasta': config[dataset]['analysis_hosts'][host],
				'virus'     : virus,
				'virus_fasta': config[dataset]['analysis_viruses'][virus],
				'analysis_condition': condition,
				'detection_mode' : detection_mode,
				'flank_region_size' : flank_region_size,
				'sensitivity_level' : sensitivity_level,
				'min_contig_length' : min_contig_length,
				'blastn_evalue_thrd': blastn_evalue_thrd,
				'similarity_thrd'   : similarity_thrd,
				'chop_read_length'  : chop_read_length,
				'minIdentity'       : minIdentity,
				'trim'              : trim,
				'merge'             : 0,
				'dedup'             : 0,
				'tool'			 : 'verse',	
			})
		i += 1
	return add_read_info(config, dataset, rows)
	

def make_vifi_rows(config, dataset):
	#### parameters for vifi ####
	i = 0
	rows = []
	# make sure the required information about the host genome has been provided
	host_file_keys = ('mappability', 'mappability_exclude', 'genes', 'exons', 'oncogenes', 'centromeres',
												'conserved_regions', 'segdup')
	hosts_to_use = []
	for host in config[dataset]['analysis_hosts'].keys():
		# check that this host is in host_info
		if host not in config[dataset]['vifi_params']['host_info']:
			print(f"host_info not provided: skipping ViFi for {host}")
			continue
		# check that all necessary files have been specified
		if not all([key in config[dataset]['vifi_params']['host_info'][host] for key in host_file_keys]):
			print(f"one or more of the required files ({host_file_keys}) for host {host} is not specfied: skipping ViFi for host {host}")
			continue
		hosts_to_use.append(host)
		
	trim_list = get_list_with_default(config[dataset]['vifi_params'], 'trim', vifi_default_trim, 'vifi_params')								

	for host, virus, trim in itertools.product(hosts_to_use, config[dataset]['analysis_viruses'].keys(), trim_list):
		condition = f"vifi{i}"
		rows.append({
				'experiment' : dataset,
				'host' 			 : host,
				'host_fasta' : config[dataset]['analysis_hosts'][host],
				'host_mappability' : config[dataset]['vifi_params']['host_info'][host]['mappability'],
				'host_mappability_exclude' : config[dataset]['vifi_params']['host_info'][host]['mappability_exclude'],				
				'host_genes' : config[dataset]['vifi_params']['host_info'][host]['genes'],				
				'host_exons' : config[dataset]['vifi_params']['host_info'][host]['exons'],				
				'host_oncogenes' : config[dataset]['vifi_params']['host_info'][host]['oncogenes'],	
				'host_centromeres' : config[dataset]['vifi_params']['host_info'][host]['centromeres'],	
				'host_conserved_regions' : config[dataset]['vifi_params']['host_info'][host]['conserved_regions'],		
				'host_segdup' : config[dataset]['vifi_params']['host_info'][host]['segdup'],																
				'virus'      : virus,		
				'virus_fasta': config[dataset]['analysis_viruses'][virus],	
				'analysis_condition': condition,
				'tool'			 : 'vifi',	
				'trim'			 : trim,
				'merge'			 : 0,
				'dedup'      : 0,	
				})
		i += 1
	return add_read_info(config, dataset, rows)


def make_polyidus_rows(config, dataset):
	#### parameters for polyidus ####
	rows = []
	i = 0
	
	trim_list = get_list_with_default(config[dataset]['polyidus_params'], 
																		'trim', 
																		polyidus_default_trim,
																		'polyidus_params')
	
	# are we trying multiple aligners?			
	aligners = get_list_with_default(config[dataset]['polyidus_params'], 
																		'aligner', 
																		polyidus_default_aligner,
																		'polyidus_params')
				
	for host, virus, aligner, trim in itertools.product(
																		config[dataset]['analysis_hosts'].keys(),
																		config[dataset]['analysis_viruses'].keys(),
																		aligners,
																		trim_list):
		# give this analysis condition a name
		condition = f"polyidus{i}"

		rows.append({
					'experiment' : dataset,
					'host' 			 : host,
					'host_fasta' : config[dataset]['analysis_hosts'][host],
					'virus'      : virus,
					'virus_fasta': config[dataset]['analysis_viruses'][virus],
					'analysis_condition': condition,
					'aligner'		 : aligner,
					'merge'			 : 0,
					'trim'       : trim,
					'dedup'      : 0,
					'tool'			 : 'polyidus',
					})	
		i += 1
	return add_read_info(config, dataset, rows)
	
def make_seeksv_rows(config, dataset):
	#### parameters for polyidus ####
	rows = []
	i = 0
	
	trim = get_list_with_default(config[dataset]['seeksv_params'], 'trim', seeksv_default_trim, 'seeksv_params')
	dedup = get_list_with_default(config[dataset]['seeksv_params'], 'dedup', seeksv_default_dedup, 'seeksv_params')
	
	trim_list = [1 if entry is True else 0 for entry in trim]
	dedup_list = [1 if entry is True else 0 for entry in dedup]	
	
	for host, virus, trim, dedup in itertools.product(
																		config[dataset]['analysis_hosts'].keys(),
																		config[dataset]['analysis_viruses'].keys(),
																		trim_list,
																		dedup_list):
		
		# give this analysis condition a name
		condition = f"seeksv{i}"

		rows.append({
					'experiment' : dataset,
					'host' 			 : host,
					'host_fasta' : config[dataset]['analysis_hosts'][host],
					'virus'      : virus,
					'virus_fasta': config[dataset]['analysis_viruses'][virus],
					'analysis_condition': condition,
					'trim'			 : trim,
					'merge'			 : 0,
					'dedup'      : dedup,
					'tool'			 : 'seeksv',
					})	
		i += 1
	return add_read_info(config, dataset, rows)


def make_vseq_rows(config, dataset):
	rows = []
	i = 0
	
	qual = get_list_with_default(config[dataset]['vseq_toolkit_params'], 'qua', vseq_default_qua, 'vseq_toolkit_params')
	lenPer = get_list_with_default(config[dataset]['vseq_toolkit_params'], 'lenPer', vseq_default_lenPer, 'vseq_toolkit_params')
	mode = get_list_with_default(config[dataset]['vseq_toolkit_params'], 'mode', vseq_default_mode, 'vseq_toolkit_params')
	vecVecFusion = get_list_with_default(config[dataset]['vseq_toolkit_params'], 'vecVecFusion', vseq_default_vecVecFusion, 'vseq_toolkit_params')
	stringencyVec = get_list_with_default(config[dataset]['vseq_toolkit_params'], 'stringencyVec', vseq_default_stringencyVec, 'vseq_toolkit_params')
	UMthresholdVec = get_list_with_default(config[dataset]['vseq_toolkit_params'], 'UMthresholdVec', vseq_default_UMthresholdVec, 'vseq_toolkit_params')
	minMapSpanVec = get_list_with_default(config[dataset]['vseq_toolkit_params'], 'minMapSpanVec', vseq_default_minMapSpanVec, 'vseq_toolkit_params')
	distVecVec = get_list_with_default(config[dataset]['vseq_toolkit_params'], 'distVecVec', vseq_default_distVecVec, 'vseq_toolkit_params')
	opVecVec = get_list_with_default(config[dataset]['vseq_toolkit_params'], 'opVecVec', vseq_default_opVecVec, 'vseq_toolkit_params')
	idenVecVec = get_list_with_default(config[dataset]['vseq_toolkit_params'], 'idenVecVec', vseq_default_idenVecVec, 'vseq_toolkit_params')
	stringencyVecGen = get_list_with_default(config[dataset]['vseq_toolkit_params'], 'stringencyVecGen', vseq_default_stringencyVecGen, 'vseq_toolkit_params')
	UMthresholdVecGen = get_list_with_default(config[dataset]['vseq_toolkit_params'], 'UMthresholdVecGen', vseq_default_UMthresholdVecGen, 'vseq_toolkit_params')
	minMapSpanVecGen = get_list_with_default(config[dataset]['vseq_toolkit_params'], 'minMapSpanVecGen', vseq_default_minMapSpanVecGen, 'vseq_toolkit_params')
	distVecGen = get_list_with_default(config[dataset]['vseq_toolkit_params'], 'distVecGen', vseq_default_distVecGen, 'vseq_toolkit_params')
	opVecGen = get_list_with_default(config[dataset]['vseq_toolkit_params'], 'opVecGen', vseq_default_opVecGen, 'vseq_toolkit_params')
	idenVecGen = get_list_with_default(config[dataset]['vseq_toolkit_params'], 'idenVecGen', vseq_default_idenVecGen, 'vseq_toolkit_params')
	clusterRange = get_list_with_default(config[dataset]['vseq_toolkit_params'], 'clusterRange', vseq_default_clusterRange, 'vseq_toolkit_params')

	# for each host, get the 'annoTable'
	annoTable = {}
	for host in config[dataset]['analysis_hosts'].keys():
		if host not in config[dataset]['vseq_toolkit_params']['host_table']:
			print(f"host_info not provided: skipping VSeq-Toolkit for {host}")
			continue	
		if host not in annoTable:
			annoTable[host] = []
		annoTable[host].append(config[dataset]['vseq_toolkit_params']['host_table'][host])				
		
	for (qual_i, lenPer_i, mode_i, vecVecFusion_i, stringencyVec_i, UMthresholdVec_i,
			minMapSpanVec_i, distVecVec_i, opVecVec_i, idenVecVec_i, stringencyVecGen_i,
			UMthresholdVecGen_i, minMapSpanVecGen_i, distVecGen_i, opVecGen_i,
			idenVecGen_i, clusterRange_i, host, virus
			) in itertools.product(
			qual, lenPer, mode, vecVecFusion, stringencyVec, UMthresholdVec,
			minMapSpanVec, distVecVec, opVecVec, idenVecVec, stringencyVecGen,
			UMthresholdVecGen, minMapSpanVecGen, distVecGen, opVecGen,
			idenVecGen, clusterRange, config[dataset]['analysis_hosts'].keys(), 
			config[dataset]['analysis_viruses'].keys()
			):
		for table in annoTable[host]:
			condition = f"vseq-toolkit{i}"
			rows.append({
				'experiment' 		: dataset,
				'host' 			 	: host,
				'host_fasta' 		: config[dataset]['analysis_hosts'][host],
				'host_table' 		: table,																
				'virus'      		: virus,		
				'virus_fasta'		: config[dataset]['analysis_viruses'][virus],	
				'analysis_condition': condition,
				'tool'				: 'vseq_toolkit',	
				'trim'				: 1,
				'merge'				: 0,
				'dedup'      		: 0,
				'qual' 				: qual_i,
				'lenPer'			: lenPer_i,
				'mode'				: mode_i,
				'vecVecFusion'		: vecVecFusion_i,
				'stringencyVec'		: stringencyVec_i,
				'UMthresholdVec'	: UMthresholdVec_i,
				'minMapSpanVec'		: minMapSpanVec_i,
				'distVecVec'		: distVecVec_i,
				'opVecVec' 			: opVecVec_i,
				'idenVecVec'		: idenVecVec_i,
				'stringencyVecGen'	: stringencyVecGen_i,
				'UMthresholdVecGen'	: UMthresholdVecGen_i,
				'minMapSpanVecGen'	: minMapSpanVecGen_i,
				'distVecGen'		: distVecGen_i,
				'opVecGen'			: opVecGen_i,
				'idenVecGen'		: idenVecGen_i,
				'clusterRange'		: clusterRange_i,
				})
			i += 1
	return add_read_info(config, dataset, rows)

def get_bool_value_from_config(config, dataset, key, default):
	if key not in config[dataset]:
		return default
		
	if config[dataset][key] is True:
		return 1
	elif config[dataset][key] is False:
		return 0
	else:
		raise ValueError(f"Boolean value for {key} in dataset {dataset} is neither True nor False.  Please specify one or the other")
		
def add_read_info(config, dataset, rows):
		# add 'read_folder', 'R1_suffix', 'R2_suffix', 'outdir', 'bam_suffix', if they are defined
 		read_folder = config[dataset].get('read_directory')
 		outdir = config[dataset].get('out_directory')
 		R1_suffix = config[dataset].get('R1_suffix')
 		R2_suffix = config[dataset].get('R2_suffix')
 		bam_suffix = config[dataset].get('bam_suffix')
 		adapter_1 = config[dataset].get('read1-adapt')
 		adapter_2 = config[dataset].get('read2-adapt')
 		
 		for row in rows:
 			row['read_folder'] = path.normpath(read_folder)
 			row['outdir'] = path.normpath(outdir)
 			row['R1_suffix'] = R1_suffix
 			row['R2_suffix'] = R2_suffix
 			row['bam_suffix'] = bam_suffix
 			row['adapter_1'] = adapter_1
 			row['adapter_2'] = adapter_2
 			
 		return rows
	
def get_samples(config):
	samples = {}
	## return a dict with datasets as key and sample names as values
	for dataset in config:
		read_dir = path.normpath(config[dataset]['read_directory'])
		samps = glob(path.join(read_dir, f"*{config[dataset]['R1_suffix']}"))
		# strip directory and suffix
		samps = [path.basename(samp)[:-len(config[dataset]['R1_suffix'])] for samp in samps]
		
		samples[dataset] = samps
		
	return samples	
		
def get_list_with_default(parent_dict, list_key, default, name):
	if list_key not in parent_dict:
		print(f"paramter {list_key} not specified for {name}: using default {default}")
		return default
		
	assert hasattr(parent_dict[list_key], '__iter__')
	return parent_dict[list_key]

