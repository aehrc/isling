#!/usr/bin/env python3

import argparse
import sys
import csv

def main(args):
	parser = argparse.ArgumentParser(description='Apply threshold to distances to calculate TP, FP, FN')
	parser.add_argument('--sim-to-found', help='Result of closest with a as simulated and b as found', required=True)
	parser.add_argument('--found-to-sim', help='Result of closest with a as found and b as sim', required=True)	
	parser.add_argument('--threshold', help='Theshold for true positive', required=True, type=int)
	parser.add_argument('--threshold-col', help='Column to apply theshold to', required=True)	
	parser.add_argument('--output', '-o', help='Output file', default='ints_summary.tsv')
	parser.add_argument('--tool', help='Tool used for analysis', required=True)
	args = parser.parse_args(args[1:])
	
	# count true positives, false positives, false negatives
	tp, fn = count_threshold(args.sim_to_found, args.threshold_col, args.threshold)
	_, fp =  count_threshold(args.found_to_sim, args.threshold_col, args.threshold)
	
	# write output
	with open(args.output, "w") as output_file:
		header = ('sim_to_found', 'found_to_sim', 'tool', 'window', 
					'coords_score_type', 'tp', 'fp', 'fn')
		results = (args.sim_to_found, args.found_to_sim, args.tool, str(args.threshold),
					args.threshold_col, str(tp), str(fp), str(fn))
					
		
		output_file.write("\t".join(header) + "\n")
		output_file.write("\t".join(results) + "\n")
		
	print(f"counted {tp} true positives, {fp} false positives, {fn} false negatives")
	print(f"saved summary to {args.output}")
	

def count_threshold(file, threshold_col, threshold):
	"""
	count number of lines with distances below and above theshold
	"""
	pos, neg = 0, 0
	tied = []
	last = (0, 0, 0, 0, 0, 0)
	with open(file, newline='') as csvfile:
		reader = csv.DictReader(csvfile, delimiter="\t")
		
		for row in reader:

			dist = float(row[f"d_{threshold_col}"])
			
			if dist < 0:
				neg += 1
			elif dist > threshold:
				neg += 1
			else:
				pos += 1
				
	return pos, neg

	

if __name__ == "__main__":
	main(sys.argv)
