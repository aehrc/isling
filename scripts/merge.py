#!/usr/bin/env python3

# merge integration sites:


# output coordinates is the range encompassing coordinates common to all ranges included in the cluster.  
# eg if a cluster consists of the host coordinates 0-1, 0-150 and 0-10, the output coordinates will be 0-1.

# file must be coordinate-sorted in host coordinates
# ie sort -k1,1 -k2,2n 
# merge first based on host coordinates, to make an intermediate cluster , then divide this cluster up into clusters based on the same process but with viral coordinates

# Two methods for clustering are available - in 'common', output clusters have some common bases in both host and virus genome
# and in 'exact', integration sites are only clustered if they have exactly the same coordinates.
# In both cases, clustered integration sites must have the same orientation (host-virus or virus-host), and 
# viral orientation (integrated in '+' or '-' orientation)


import argparse
import sys
import csv
import pdb
import pprint

def main(args):

	parser = argparse.ArgumentParser(description = "merge integration sites based on host and virus coordinates")
	parser.add_argument('--input', '-i', help = 'input file (from postprocessing)', required=True)
	parser.add_argument('--output', '-o', help = 'output file', default='merged.txt')
	parser.add_argument('--min-n', '-n', help = 'minimum number of reads to retain a cluster', default=1, type=int)
	parser.add_argument('--cluster-method', '-c', help = 'method for clustering', choices = {'common', 'exact'}, default='common')
	args = parser.parse_args()
	
	
	#open input, output files
	with open(args.input, 'r', newline = '') as infile, open(args.output, 'w', newline = '') as outfile:
		reader = csv.DictReader(infile, delimiter = '\t')
		writer = csv.writer(outfile, delimiter = '\t')

		# Write header to file
		header = ('Chr', 'IntStart', 'IntStop', 'JunctionType', 'Virus', 'VirusStart', 'VirusStop', 'VirusOri', 'nChimeric', 'nDiscordant', 'SiteID', 'ReadIDs')
		writer.writerow(header)	
		
		# get first row to initalize
		try:
			row = next(reader)

		# if file is empty, just quit
		except StopIteration:
			return
		
		row = prune_row(row)
		n_line = 1		
		# variables for checking sorting
		curr = reset_curr(row)
		seen_host_chrs = {row['Chr']}
		
		# variables for merging clusters
		n_clust = 0
		clust = reset_clust(row)
	
		for row in reader:

			row = prune_row(row)

			# check sorting
			seen_host_chrs = check_sorted(row, curr, n_line, seen_host_chrs)
			
			# set current values
			curr = reset_curr(row)

			# initial merging only based on host
			merge = check_overlap(row, clust, 'host')

			if merge:
				clust = merge_row(clust, row)
						
			else:
				
				clusts = process_initial_cluster(clust, method = args.cluster_method)
				
				# otherwise, write previous cluster to outfile
				write_clusters(writer, clusts, n_clust, args.min_n)
							
				# start new cluster			
				n_clust += len(clusts)					
				clust = reset_clust(row)		
			
			n_line += 1 
		
		# write last cluster to file
		clusts = process_initial_cluster(clust, method = args.cluster_method)
		write_clusters(writer, clusts, n_clust, args.min_n)
		
		print(f"read {n_line} reads from file {args.input}")
		print(f"grouped into {n_clust} clusters using method '{args.cluster_method}'")
		print(f"saved clusters to file {args.output}")

def process_initial_cluster(clust, method = 'common'):
	"""
	Process initial cluster (based on host coordinates) to produce further clusters.
	 
	In method 'common', reads belong in the same cluster if they have the same host and virus orientation
	and they also have common coordinates in both host and virus
	
	In 'exact', reads belong in the same cluster if they ahve the same host and virus orientation
	and they have the same coordiantes in both host and virus
	"""
	assert method in {'common', 'exact'}
	
	# if there's only one row in this cluster, we're already done
	if len(clust['ints']) == 1:
		return [clust]
	
	# create new cluster from first integration
	clusts = [reset_clust(clust['ints'].pop())]
	
	# for each integration in this cluster, check if it belongs in any of the current clusters
	for row in clust['ints']:
		merge = False
		for c in clusts:
			
			# if we already merged this row with a cluster, skip
			if merge:
				continue
		
			# check if we should merge row in to cluster c
			if method == 'common' and merge_common_coords(c, row):
				merge = True
				c = merge_row(c, row, method = 'common')
				continue
			elif method == 'exact' and merge_same_coords(c, row):
				merge = True
				c = merge_row(c, row, method = 'common') # method doesn't matter because we already know the coordinates are the same
		
		# if we didn't merge this row with a cluster, we need a new cluster
		if not merge:
			clusts.append(reset_clust(row))
			
	return clusts
			
				
	
def merge_common_coords(clust, row):
	"""
	Check if both host and virus have an overlap
	"""
	if not check_overlap(row, clust, 'host', True):
		return False
	if not check_overlap(row, clust, 'virus', True):
		return False
		
	return True
	
def merge_same_coords(clust, row):
	"""
	Check if both host and virus have same coordinates
	"""	
	if not check_same_coordinates(row, clust, 'host', True):
		return False
	if not check_same_coordinates(row, clust, 'virus', True):
		return False
	
	return True



def get_overlap(start1, stop1, start2, stop2):
	"""
	return the coordinates of the overlaped bases between two intervals.
	eg if one integration has coordinates (50, 60) and another has coordinates (40, 100), the comon/overlapped 
	coordinates would be (50, 60)
	"""
	assert overlap('_', start1, stop1, '_', start2, stop2)
	
	# get the higher of the start coordinates
	if start1 > start2:
		start = start1
	else:
		start = start2
	# get the lower of the stop coordinates
	if stop1 < stop2:
		stop = stop1
	else:
		stop = stop2
		
	return (start, stop)

	
def reset_curr(row):
		"""
		reset curr dict with data from input row
		"""
		return {'host_chr'		: row['Chr'],
						'host_start'	: int(row['IntStart']),
						'virus'				: row['VirusRef'],
						'virus_start'	: int(row['VirusStart'])
					}
			
def reset_clust(row):
	"""
	reset clust dict with data from input row
	"""			
	clust = {
		'host_chr'		: row['Chr'],
		'host_coords'	: [int(row['IntStart']), int(row['IntStop'])], # stores range encompassing all host coordinates in cluster
		'oris'				: {row['Orientation']},
		'virus'				: row['VirusRef'],
		'virus_coords': [int(row['VirusStart']), int(row['VirusStop'])],# stores range encompassing all virus coordinates in cluster
		'virus_ori'		: {row['VirusOrientation']},
		'reads'				: {row['ReadID']},
		'n_chimeric'  : 1 if row['OverlapType'] != 'discordant' else 0,
		'n_discord'		: 1 if row['OverlapType'] == 'discordant' else 0,
		'ints': [row] # store all integrations currently in cluster
	}

	return clust

def write_clusters(outcsv, clusts, clust_id, min_n):
	"""
	write cluster row to outfile
	"""
	for clust in clusts:
		if len(clust['reads']) >= min_n:
			row = (clust['host_chr'], clust['host_coords'][0], clust['host_coords'][1], clust['oris'].pop(),
						clust['virus'], clust['virus_coords'][0], clust['virus_coords'][1],  clust['virus_ori'].pop(),
						clust['n_chimeric'], clust['n_discord'],
						clust_id, ",".join(clust['reads']))
			outcsv.writerow(row)
			clust_id += 1

def check_overlap(row, clust, ref, check_ori = False):
	"""
	check if there's overlap in both host intervals defined in row, and intervals defined in clust
	'ref' is either 'host' or 'virus', and indicates whether to check for overlap in host or virus
	"""
	assert ref in {'host', 'virus'}
	assert check_ori is True or check_ori is False
	if check_ori and ref == 'host':
		assert len(clust['oris']) == 1
	elif check_ori and ref == 'virus':
		assert len(clust['virus_ori']) == 1
	
	if ref == 'host':
		if check_ori and (row['Orientation'] not in clust['oris']):
			return False
		return overlap(row['Chr'], int(row['IntStart']), int(row['IntStop']), 
									clust['host_chr'], clust['host_coords'][0], clust['host_coords'][1])
	else:
		if check_ori and (row['VirusOrientation'] not in clust['virus_ori']):
			return False		
		return overlap(row['VirusRef'], int(row['VirusStart']), int(row['VirusStop']), 
									clust['virus'], clust['virus_coords'][0], clust['virus_coords'][1])

def check_same_coordinates(row, clust, ref, check_ori = True):
	"""
	check if both host intervals defined in row, and intervals defined in clust have same start and stop coordinates
	'ref' is either 'host' or 'virus', and indicates whether to check for overlap in host or virus
	"""
	assert ref in {'host', 'virus'}
	assert check_ori is True or check_ori is False
	if check_ori and ref == 'host':
		assert len(clust['oris']) == 1
	elif check_ori and ref == 'virus':
		assert len(clust['virus_ori']) == 1
		
	if ref == 'host':
		if check_ori and (row['Orientation'] not in clust['oris']):
			return False
		if row['Chr'] != clust['host_chr']:
			return False
		if int(row['IntStart']) != clust['host_coords'][0]:
			return False
		if int(row['IntStop']) != clust['host_coords'][1]:
			return False
	else:
		if check_ori and (row['VirusOrientation'] not in clust['virus_ori']):
			return False
		if row['VirusRef'] != clust['virus']:
			return False
		if int(row['VirusStart']) != clust['virus_coords'][0]:
			return False
		if int(row['VirusStop']) != clust['virus_coords'][1]:
			return False
			
	return True
	
def check_sorted(row, curr, n_line, seen_host_chrs):
	"""
	check that position in host is the same as previous, or increasing
	"""
	
	# check coordinates are valid
	assert int(row['IntStart']) <= int(row['IntStop'])
	assert int(row['VirusStart']) <= int(row['VirusStop'])
	
	# check sorting
	# same chromosome as current
	if row['Chr'] == curr['host_chr']:
		# check for position before current position
		if int(row['IntStart']) < curr['host_start']:
			raise ValueError(f'Input file must be sorted (ie `sort -k1,1 -k2,2n`): problem with row {n_line}')				
									
	# chromosome is different	
	else:
		# if chomosome is different but we've seen it before
		if row['Chr'] in seen_host_chrs:
			raise ValueError(f'Input file must be sorted (ie `sort -k1,1 -k2,2n`): problem with row {n_line}')
		else:
			# add to list of seen chromocomes
			seen_host_chrs.add(row['Chr'])
					
	return seen_host_chrs	
	
def overlap(ref1, start1, stop1, ref2, start2, stop2):
	"""
	check if intervals ref1:start1-stop1 and ref2:start2-stop2 overlap
	"""
	assert stop1 >= start1
	assert stop2 >= start2
	
	# different chromosome
	if ref1 != ref2:
		return False
	# interval 1 completely to the left of interval 2
	if start1 > stop2:
		return False
	# interval 2 completely to the left of interval 1
	if start2 > stop1:
		return False
	return True
	
def merge_row(clust, row, method = 'all_bases'):
	"""
	merge integration defined in row with current cluster
	
	two methods for merging coordinates:
		'all_bases': coordinates of cluster encompass full range of coordinates occupied by both cluster and row
		'common': coordinates of cluster are only those bases common to the cluster and the row
	"""
	
	assert method in {'all_bases', 'common'}
	assert clust['host_chr'] == row['Chr']
	
	clust['reads'].add(row['ReadID'])
	clust['ints'].append(row)
	clust['oris'].add(row['Orientation'])
	clust['virus_ori'].add(row['VirusOrientation'])
	
	if row['OverlapType'] == 'discordant':
		clust['n_discord'] += 1
	else:
		clust['n_chimeric'] += 1
					
	# re-calculate cluster coordinates
	if method == 'all_bases':
		clust['host_coords'][0] = min(int(row['IntStart']), clust['host_coords'][0])
		clust['host_coords'][1] = max(int(row['IntStop']), clust['host_coords'][1])
		clust['virus_coords'][0] = min(int(row['VirusStart']), clust['virus_coords'][0])
		clust['virus_coords'][1] = max(int(row['VirusStop']), clust['virus_coords'][1])
	else:
		clust['host_coords'] = get_overlap(int(row['IntStart']), int(row['IntStop']), 
																				clust['host_coords'][0],  clust['host_coords'][1])
		clust['virus_coords'] = get_overlap(int(row['VirusStart']), int(row['VirusStop']), 
																				clust['virus_coords'][0],  clust['virus_coords'][1])
		
		
	return clust

def prune_row(row):
	
	# remove extraneous information
	tmp = {}
	for key in row:

		if key in {'Chr', 'IntStart', 'IntStop', 'Orientation', 'VirusRef', 'VirusStart', 'VirusStop', 'VirusOrientation', 'OverlapType', 'ReadID'}:
			tmp[key] = row[key]			
	return tmp
	
	
if __name__ == "__main__":
	main(sys.argv)

