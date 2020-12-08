#!/usr/bin/env python3

# merge integration sites in various ways:

# --type exact: merge only sites with exactly the same coordinates in both host and virus
# --type overlap: merge any sites with overlapping coordinates in both host and virus,  output the range of coordinates
# --type distance: merge any sitest within specified of another site in both host and virus
# if type is 'distance', must also specify -d <merge_distance>

# if type is 'overlap' or 'distance', output coordinates is the range encompassing coordinates common to all ranges included in the cluster.  eg if a cluster consists of the host coordinates 0-1, 0-150 and 0-10, the output coordinates will be 0-1.
# if there are no coordinates common to all included ranges, (STILL NEED TO SORT THIS OUT!!!)

# file must be coordinate-sorted in host coordinates
# ie sort -k1,1 -k2,2n 
# merge first based on host coordinates, to make an intermediate cluster , then divide this cluster up into clusters based on the same process but with viral coordinates

# for now, don't consider orientation of integration site (host/virus or virus/host), but just merge based on coordinates


import argparse
import sys
import csv
import pdb
import pprint

def main(args):

	parser = argparse.ArgumentParser(description = "merge integration sites based on host and virus coordinates")
	parser.add_argument('--input', '-i', help = 'input file (from postprocessing)', required=True)
	parser.add_argument('--output', '-o', help = 'output file', default='merged.txt')
	parser.add_argument('--distance', '-d', help = 'distance for merging', type=int, default = 0)
	parser.add_argument('--min-n', '-n', help = 'minimum number of reads to retain a cluster', default=1, type=int)
	args = parser.parse_args()
	
	
	#open input, output files
	with open(args.input, 'r', newline = '') as infile, open(args.output, 'w', newline = '') as outfile:
		reader = csv.DictReader(infile, delimiter = '\t')
		writer = csv.writer(outfile, delimiter = '\t')

		# Write header to file
		header = ('Chr', 'IntStart', 'IntStop', 'Virus', 'VirusStart', 'VirusStop', 'nChimeric', 'nDiscordant', 'SiteID', 'ReadIDs')
		writer.writerow(header)	

		# Check if input file is empty. If empty write file only with header

		# get first row to initialize 
		try:
			row = next(reader)
			row = prune_row(row)
			n_line = 1
		except StopIteration:
			return
		
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
			merge = check_overlap(row, clust, args.distance, 'host')

			if merge:
				clust = merge_row(clust, row)
						
			else:
				# process this cluster - cluster using viral coordinates
				clusts = cluster_virus(clust, args.distance)
				
				# get common coordinates for output - this might result in further splitting clusters if there are no common coordinates
				# for a particular cluster (ie if d was too large and resulted in merging of independent integrations)		
				clusts = add_common_coords(clusts)
				
				# otherwise, write previous cluster to outfile
				write_clusters(writer, clusts, n_clust, args.min_n)
							
				# start new cluster			
				n_clust += len(clusts)					
				clust = reset_clust(row)		
			
			n_line += 1 
		
		# write last cluster to file
		clusts = cluster_virus(clust, args.distance)
		clusts = add_common_coords(clusts)
		write_clusters(writer, clusts, n_clust, args.min_n)

def add_common_coords(clusts):
	"""
	add common coordinates to each cluster - the coordinates (host and virus) that are common to all integrations within the cluster
	eg if one integration has coordinates (50, 60) and another has coordinates (40, 100), the comon coordinates would be (50, 60)
	
	if there are no common coordinates, try to group into multiple clusters that each have some common coordinates
	"""	
	done = []
	
	# try each cluster in turn
	for clust in clusts:	
	
		# if there's only one integration in this cluster
		if len(clust['ints']) == 1:
			done.append(clust)
			continue

		# sort integrations by host coords
		clust['ints'].sort(key = lambda row: int(row['IntStart']))
		
		# keep getting common coordinates in host until we find a non-overlap
		breakpoint = -1
		for i, row in enumerate(clust['ints']):
			try:
				clust['host_coords'] = get_overlap(clust['host_coords'][0], clust['host_coords'][1], 
																						int(row['IntStart']), int(row['IntStop']))
			
			# if non-overlap, we need to split the cluster in two
			except AssertionError:
				breakpoint = i
				break
				
		if breakpoint != -1	:
			# retry these two new clusters
			done += add_common_coords(split_clust(clust, breakpoint))
			continue
				
		# sort by virus coords
		clust['ints'].sort(key = lambda row: int(row['VirusStart']))		
		
		# keep getting common coordinates in virus until we find a non-overlap
		for i, row in enumerate(clust['ints']):
			try:
				clust['virus_coords'] = get_overlap(clust['virus_coords'][0], clust['virus_coords'][1], 
																						int(row['VirusStart']), int(row['VirusStop']))
			# if non-overlap, we need to split the cluster in two
			except AssertionError:
				breakpoint = i
				break
			
		if breakpoint != -1	:
			# retry these two new clusters
			done +=  add_common_coords(split_clust(clust, breakpoint))
			continue		

		# if we got all the way through, add this cluster
		done.append(clust)
		
	return done

def split_clust(clust, i):
	"""
	split a cluster into two at row number i
	"""
	clust1 = reset_clust(clust['ints'][0])
	for j in range(i):
			clust1 = merge_row(clust1, clust['ints'][j])
	
	clust2 = reset_clust(clust['ints'][i])
	for j in range(i, len(clust['ints'])):
		clust2 = merge_row(clust2, clust['ints'][j])
			
	return [clust1, clust2]


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

def cluster_virus(clust, d):
	"""
	process initial cluster (based on host coordinates) - divide into smaller clusters based on distance in virus
	ideally we should sort integrations, but assume there won't be too many in each cluster so it's probably not necessary
	"""
	# if there's only one row in this cluster, we're already done
	if len(clust['ints']) == 1:
		return [clust]	
	
	# otherwise, sort the integrations	
	clust['ints'].sort(key = lambda row: int(row['VirusStart']))
	
	# start with first integration in cluster
	clusts = [reset_clust(clust['ints'][0])]
	
	# check remaining rows
	curr_clust = 0
	for row in clust['ints'][1:]:
		
		# check for overlap with current cluster
		if check_overlap(row, clusts[curr_clust], d, 'virus') and check_overlap(row, clusts[curr_clust], d, 'host'):
			clusts[curr_clust] = merge_row(clusts[curr_clust], row)
		
		# if no overlap, start new cluster	
		else:
			clusts.append(reset_clust(row))
			curr_clust += 1
		
			
	
	return clusts
	
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
		'virus'				: row['VirusRef'],
		'virus_coords'	: [int(row['VirusStart']), int(row['VirusStop'])],# stores range encompassing all virus coordinates in cluster
		'reads'				: [row['ReadID']],
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
			row = (clust['host_chr'], clust['host_coords'][0], clust['host_coords'][1], 
						clust['virus'], clust['virus_coords'][0], clust['virus_coords'][1], 
						clust['n_chimeric'], clust['n_discord'],
						clust_id, ",".join(clust['reads']))
			outcsv.writerow(row)
			clust_id += 1



def check_overlap(row, clust, dist, ref):
	"""
	check if there's overlap in both host intervals defined in row, and intervals defined in clust
	'ref' is either 'host' or 'virus', and indicates whether to check for overlap in host or virus
	"""
	assert ref in {'host', 'virus'}
	
	if ref == 'host':
		return overlap(row['Chr'], int(row['IntStart']) - dist, int(row['IntStop']) + dist, 
									clust['host_chr'], clust['host_coords'][0], clust['host_coords'][1])
	else:
		return overlap(row['VirusRef'], int(row['VirusStart']) - dist, int(row['VirusStop']) + dist, 
									clust['virus'], clust['virus_coords'][0], clust['virus_coords'][1])
	
def check_sorted(row, curr, n_line, seen_host_chrs):
	"""
	check that position in host is the same as previous, or increasing
	"""
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
	
def merge_row(clust, row):
	"""
	merge integration defined in row with current
		'host_chr'		: row['Chr'],
		'host_coords'	: [int(row['IntStart']), int(row['IntStop'])], # stores range encompassing all host coordinates in cluster
		'virus'				: row['VirusRef'],
		'virus_coords'	: [int(row['VirusStart']), int(row['VirusStop'])],# stores range encompassing all virus coordinates in cluster
		'reads'				: [row['ReadID']],
		'n_chimeric'  : 1 if row['OverlapType'] != 'discordant' else 0,
		'n_discord'		: 1 if row['OverlapType'] == 'discordant' else 0,
		'ints': [row] # store all integrations currently in cluster
	"""
	
	clust['reads'].append(row['ReadID'])
	clust['ints'].append(row)
	
	if row['OverlapType'] == 'discordant':
		clust['n_discord'] += 1
	else:
		clust['n_chimeric'] += 1
					
	# extend cluster start and stop, if necessary
	clust['host_coords'][0] = min(int(row['IntStart']), clust['host_coords'][0])
	clust['host_coords'][1] = max(int(row['IntStop']), clust['host_coords'][1])
	clust['virus_coords'][0] = min(int(row['VirusStart']), clust['virus_coords'][0])
	clust['virus_coords'][1] = max(int(row['VirusStop']), clust['virus_coords'][1])
		
	return clust

def prune_row(row):
	# remove extraneous information
	tmp = {}
	for key in row:
		if key in {'Chr', 'IntStart', 'IntStop', 'VirusRef', 'VirusStart', 'VirusStop', 'OverlapType', 'ReadID'}:
			tmp[key] = row[key]
					
	return tmp
	
	
if __name__ == "__main__":
	main(sys.argv)

