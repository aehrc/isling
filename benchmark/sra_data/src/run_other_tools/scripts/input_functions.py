import pdb
import os

def analysis_df_value(wildcards, analysis_df, column_name):
	"""Get resources required for a dataset (but not tool-specific)"""
	#return analysis_df.loc[(analysis_df['experiment'] == wildcards.dset).idxmax(), column_name] 
	row = analysis_df[(analysis_df['experiment'] == wildcards.dset)].index[0]
	return analysis_df.loc[row, column_name]

def analysis_df_tool_value(wildcards, analysis_df, tool, column_name):
	"""Get resources required for a dataset and tool"""
	row = analysis_df[(analysis_df['experiment'] == wildcards.dset) & (analysis_df['tool'] == tool)].index[0]
	return analysis_df.loc[row, column_name]



def get_input_reads(wildcards, analysis_df, read_num):
	assert read_num in (1, 2)
	if read_num == 1:
		return f"{analysis_df_value(wildcards, analysis_df, 'read_folder')}/{wildcards.samp}{analysis_df_value(wildcards, analysis_df,  'R1_suffix')}"
	if read_num == 2:
		return f"{analysis_df_value(wildcards, analysis_df,  'read_folder')}/{wildcards.samp}{analysis_df_value(wildcards, analysis_df,  'R2_suffix')}"	

def get_vifi_resource(wildcards, analysis_df, resource_name):
	"""Get resources required for vifi"""
	host_idx = analysis_df[(analysis_df['host'] == wildcards.host) & (analysis_df['tool'] == 'vifi')].index[0]
	return analysis_df.loc[host_idx, resource_name]
	
def get_reads(wildcards, analysis_df, rules, read_num):
	assert read_num in (1, 2)
	if analysis_df_value(wildcards, analysis_df, 'trim') == 1:
		if read_num == 1:
			return rules.trim.output.proc_r1
		else:
			return rules.trim.output.proc_r2
	else:
		if read_num == 1:
			return get_input_reads(wildcards, analysis_df, 1)
		else:
			return get_input_reads(wildcards, analysis_df, 2)
			

def get_vifi_resource(wildcards, analysis_df, resource_name):
	"""Get resources required for vifi"""
	host_idx = analysis_df[(analysis_df['host'] == wildcards.host) & (analysis_df['tool'] == 'vifi')].index[0]
	return analysis_df.loc[host_idx, resource_name]
	


def resources_list_with_min_and_max(file_name_list, attempt, mult_factor=2, minimum = 100, maximum = 50000):

	resource = int(sum([os.stat(file).st_size/1e6 for file in file_name_list]) * mult_factor * attempt) 
	
	resource = min(maximum, resource)
	
	return max(minimum, resource)


def tool_datsets_regex(wildcards, analysis_df):
	# get tool from wildcards
	tool = analysis_df_value(wildcards, analysis_df, 'tool')

	# get all analysis conditions for this tool
	return "|".join(analysis_df.loc[analysis_df['tool'] == tool, 'analysis_condition'])
	
	

