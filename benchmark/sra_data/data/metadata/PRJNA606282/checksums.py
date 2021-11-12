#!/usr/bin/env python3

#usage: python3 checksums.py <ENA_file>

import sys
import os

ena=sys.argv[1]

with open(ena, 'r') as enafile:
	for line in enafile:
		parts = line.strip().split('\t')
		sums = parts[6].split(";")
		names = parts[7].split(";")
		
		for i in range(len(sums)):
			print(f"{sums[i]}\t{os.path.basename(names[i])}")
