#!/usr/bin/env python3

# merge integration sites in various ways:

# --type exact: merge only sites with exactly the same coordinates in both host and virus
# --type overlap: merge any sites with overlapping coordinates in both host and virus,  output the range of coordinates
# --type distance: merge any sitest within specified of another site in both host and virus
# if type is 'distance', must also specify -d <merge_distance>

# if type is 'overlap' or 'distance', output coordinates is the range encompassing coordinates of all chimeric reads in the pair.  if there are no chimeric reads, then the range encompassing all coordinates will be output

# file must be coordinate-sorted in both host and virus coordinates
# ie sort -k1,1 -k2,2n -k4,4 -k5,5n

# for now, don't consider orientation of integration site (host/virus or virus/host), but just merge based on coordinates


import argparse
import sys
import csv
import pdb

def main(args):

	parser = argparse.ArgumentParser(description = "merge integration sites based on host and virus coordinates")
	parser.add_argument('--input', '-i', help = 'input file (from postprocessing)', required=True)
	parser.add_argument('--output', '-o', help = 'output file', default='merged.txt')
	parser.add_argument('--type', '-t', help = 'type of merging to do', choices = ("exact", "overlap", "distance"), required=True)
	parser.add_argument('--distance', '-d', help = 'distance for merging', type=int)
	parser.add_argument('--min-n', '-n', help = 'minimum number of reads to retain a cluster', default=1, type=int)
	args = parser.parse_args()
	
	# if type == 'distance', check -d was specified
	if args.type == 'distance':
		if args.distance is None:
			raise ValueError("Please specify positive integer -d when type is 'distance'")
		if args.distance < 0:
			raise ValueError("Please specify positive integer -d when type is 'distance'")
	
	# type 'overlap' is equivalent to type 'distance' with -d 0
	if args.type == 'overlap':
		args.distance = 0
		args.type = 'distance'
	
	#open input, output files
	with open(args.input, 'r', newline = '') as infile, open(args.output, 'w', newline = '') as outfile:
		reader = csv.DictReader(infile, delimiter = '\t')
		writer = csv.writer(outfile, delimiter = '\t')

		# Write header to file
		header = ('Chr', 'IntStart', 'IntStop', 'Virus', 'VirusStart', 'VirusStop', 'nChimeric', 'nDiscordant', 'SiteID', 'ReadIDs')
		writer.writerow(header)	

		# Check if input file is empty. If empty write file only with header
		is_empty = False
		for row in reader:
			is_empty = True
		if not is_empty:
			return None

		# get first row to initialize 
		row = next(reader)
		n_line = 1
		
		# variables for checking sorting
		curr = reset_curr(row)
		seen_host_chrs = {row['Chr']}
		seen_viruses = {row['VirusRef']}
		
		# variables for merging clusters
		n_clust = 0
		clust = reset_clust(row)
	
		for row in reader:
		
			# check sorting
			# same chromosome as current
			if row['Chr'] == curr['host_chr']:
				# check for position before current position
				if int(row['IntStart']) < curr['host_start']:
					raise ValueError(f'Input file must be sorted (ie `sort -k1,1 -k2,2n -k4,4 -k5,5n`): problem with row {n_line}')
				
				# position same 
				elif int(row['IntStart']) == curr['host_start']:
					seen_viruses = check_virus_sorted(row, curr, seen_viruses, n_line)
				
				# position greater than current position
				else:
					# if the position is not the same, reset the list of seen viruses
					seen_viruses = set()					
									
			# chromosome is different	
			else:
				# if chomosome is different but we've seen it before
				if row['Chr'] in seen_host_chrs:
					raise ValueError(f'Input file must be sorted (ie `sort -k1,1 -k2,2n -k4,4 -k5,5n`): problem with row {n_line}')
				else:
					# add to list of seen chromocomes
					seen_host_chrs.add(row['Chr'])
					seen_viruses = set()
					
			# set current values
			curr = reset_curr(row)
			
			if args.type == 'exact':
				# check if we should merge this row with current cluster
				merge = check_merge_exact(row, clust)
						
				# if we're merging this read
				if merge:
					clust['reads'].append(row['ReadID'])
					if row['OverlapType'] == 'discordant':
						clust['n_discord'] += 1
					else:
						clust['n_chimeric'] += 1						
				
				else:
					# otherwise, write previous cluster to outfile
					write_cluster_simple(writer, clust, n_clust, args.min_n)
							
					# start new cluster			
					n_clust += 1					
					clust = reset_clust(row)
			
			elif args.type == 'distance':
				# check if we should merge this row with current cluster
				merge = check_merge_overlap(row, clust, args.distance)

				if merge:
					clust['reads'].append(row['ReadID'])
					if row['OverlapType'] == 'discordant':
						clust['n_discord'] += 1
					else:
						clust['n_chimeric'] += 1
					
					# if this read is chimeric, but cluster previously contained only discordant pairs
					# use the coorindates of this chimeric read instead of the di
					
					# extend cluster stop, if necessary
					# don't need to extend start because input is sorted
					clust['host_stop'] = max(int(row['IntStop']), clust['host_stop'])
					clust['virus_stop'] = max(int(row['VirusStop']), clust['virus_stop'])
					
					# if this read is chimeric, also adjust chimeric start and stop if necessary
					if row['OverlapType'] != 'discordant' and not clust['contains_chimeric']:
						clust['chimeric_host_start'] = int(row['IntStart'])
						clust['chimeric_host_stop'] = int(row['IntStop'])
						clust['chimeric_virus_start'] = int(row['VirusStart'])
						clust['chimeric_virus_stop'] = int(row['VirusStop'])
						clust['contains_chimeric'] = True
					elif row['OverlapType'] != 'discordant' and clust['contains_chimeric']:
						clust['chimeric_host_stop'] = max(int(row['IntStop']), clust['chimeric_host_stop'])
						clust['chimeric_virus_stop'] = max(int(row['VirusStop']), clust['chimeric_virus_stop'])				
						
				else:
					# otherwise, write previous cluster to outfile
					write_cluster(writer, clust, n_clust, args.min_n)
							
					# start new cluster			
					n_clust += 1					
					clust = reset_clust(row)		
			
			n_line +=1 
		
		# write last cluster to file
		if args.type == 'exact':
			write_cluster_simple(writer, clust, n_clust, args.min_n)
		else:
			write_cluster(writer, clust, n_clust, args.min_n)
			
	
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
		'host_start'	: int(row['IntStart']),
		'host_stop'		: int(row['IntStop']),
		'virus'				: row['VirusRef'],
		'virus_start'	: int(row['VirusStart']),
		'virus_stop'	: int(row['VirusStop']),
		'reads'				: [row['ReadID']],
		'contains_chimeric': (row['OverlapType'] != 'discordant'),
		'n_chimeric'  : 1 if row['OverlapType'] != 'discordant' else 0,
		'n_discord'		: 1 if row['OverlapType'] == 'discordant' else 0
	}
	if row['OverlapType'] != 'discordant':
		clust['chimeric_host_start'] = clust['host_start']
		clust['chimeric_host_stop'] = clust['host_stop']
		clust['chimeric_virus_start'] = clust['virus_start']
		clust['chimeric_virus_stop'] = clust['virus_stop']
		
	return clust


def write_cluster_simple(outcsv, clust, clust_id, min_n):
	"""
	write cluster row to outfile
	"""
	if len(clust['reads']) >= min_n:
		row = (clust['host_chr'], clust['host_start'], clust['host_stop'], 
					clust['virus'], clust['virus_start'], clust['virus_stop'], 
					clust['n_chimeric'], clust['n_discord'],
					clust_id, ",".join(clust['reads']))
		outcsv.writerow(row)

def write_cluster(outcsv, clust, clust_id, min_n):
	"""
	write cluster row to outfile
	"""
	if len(clust['reads']) >= min_n:	
		if clust['contains_chimeric']:
			row = (clust['host_chr'], clust['chimeric_host_start'], clust['chimeric_host_stop'], 
							clust['virus'], clust['chimeric_virus_start'], clust['chimeric_virus_stop'], 
							clust['n_chimeric'], clust['n_discord'],
							clust_id, ",".join(clust['reads']))
			outcsv.writerow(row)
		else:
			write_cluster_simple(outcsv, clust, clust_id, min_n)

	
def check_merge_exact(row, clust):
	"""
	check for an exact match between coordinates from this row and current cluster coordinates
	return True for exact match and False otherwise
	"""
	if row['Chr'] != clust['host_chr']:
		return False
	if int(row['IntStart']) != clust['host_start']:
		return False
	if int(row['IntStop']) != clust['host_stop']:
		return False
	if row['VirusRef'] != clust['virus']:
		return False
	if int(row['VirusStart']) != clust['virus_start']:
		return False
	if int(row['VirusStop']) != clust['virus_stop']:
		return False
	
	return True

def check_merge_overlap(row, clust, dist):
	"""
	check if there's overlap in both host and virus between intervals defined in row and intervals defined in clust
	"""
	# overlap in host
	if not overlap(row['Chr'], int(row['IntStart']) - dist, int(row['IntStop']) + dist, 
									clust['host_chr'], clust['host_start'], clust['host_stop']):	
		return False
	# overlap in virus
	if not overlap(row['VirusRef'], int(row['VirusStart']) - dist, int(row['VirusStop']) + dist, 
									clust['virus'], clust['virus_start'], clust['virus_stop']):
		return False
		
	return True

	
def check_virus_sorted(row, curr, seen_viruses, n_line):
	"""
	check that position in virus is the same as previous, or increasing
	"""
	# if the same as current virus
	if row['VirusRef'] == curr['virus']:
		# check for virus position before current position
		if int(row['VirusStart']) < curr['virus_start']:
			raise ValueError(f'Input file must be sorted (ie `sort -k1,1 -k2,2n -k4,4 -k5,5n`): problem with row {n_line}')
	# if not the same as current virus
	else:
		# check we haven't already seen this virus
		if row['VirusRef'] in seen_viruses:
			raise ValueError(f'Input file must be sorted (ie `sort -k1,1 -k2,2n -k4,4 -k5,5n`): problem with row {n_line}')
		else:
			seen_viruses.add(row['VirusRef'])			
		
	return seen_viruses
	
	
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

if __name__ == "__main__":
	main(sys.argv)

