
def make_reference_dict(config):
	# construct dictionary with reference names as keys and reference fastas as values
	host_names = { config[dataset]['host_name']:config[dataset]['host_fasta'] for dataset in config }
	virus_names = { config[dataset]['virus_name']:config[dataset]['virus_fasta'] for dataset in config }
	return {**host_names, **virus_names}
	
	
