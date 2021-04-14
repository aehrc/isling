
def make_reference_dict(toDo):
	# construct dictionary with reference names as keys and reference fastas as values
	
	ref_names = {}
	for index, row in toDo.iterrows():
		ref_names[row['host']] = row['host_fasta']
		ref_names[row['virus']] = row['virus_fasta']
	
	return ref_names
	
	
