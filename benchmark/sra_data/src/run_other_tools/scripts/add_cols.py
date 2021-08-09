#!/usr/bin/env python3

import sys

header = ("chr_1", "start_1", "stop_1", "chr_2", "start_2", "stop_2", 
			"shortest", "coords_mean", "coords_min", "midpoint")

def main(stdin):

	print("\t".join(header))
	
	for line in stdin:
		parts = line.strip().split()
		
		# if there is no closest integration
		if parts[3] == ".":
			coords_mean, coords_min, midpoint = -1, -1, -1
		else:
			# difference between each start and stop
			diff_1 = abs(int(parts[1])- int(parts[4]))
			diff_2 = abs(int(parts[2])- int(parts[5]))
			coords_mean = (diff_1 + diff_2)/2
			coords_min = min(diff_1, diff_2)
			
			# midpoint
			mid_1 = (int(parts[2]) + int(parts[1]))/2
			mid_2 = (int(parts[5]) + int(parts[4]))/2
			
			midpoint = abs(mid_2 - mid_1)
			
		line = parts + [coords_mean, coords_min, midpoint]
		line = [str(i) for i in line]
		print("\t".join(line))
			

if __name__ == "__main__":
	main(sys.stdin)
