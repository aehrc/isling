# usage: python3 combine_ucsc.py <file1> <file2> <file3> ... <out>

import sys
import csv
import re
import subprocess as sp
import pdb

# have to have at least one file to combined and an outfile name
def main():
	assert len(sys.argv) >= 3

	print(f"combining {len(sys.argv) - 2} files")
	outfile = sys.argv[-1]
	tmpfile = outfile + ".tmp"

	# combine files into one
	with open(tmpfile, 'w') as csvfile:

		out_writer = csv.writer(csvfile, delimiter='\t')
	
		for infile in sys.argv[1:-1]:
			print(f"working on {infile}")
			with open(infile, 'r', newline='') as incsv:
			
				# get track name
				try:
					track_header = next(incsv)
					next(incsv) # skip header
				except StopIteration:
					print(f"file {infile} is empty! skipping")
					continue
				
				track_name = get_track_name(track_header)
				
			
				in_reader = csv.reader(incsv, delimiter='\t')
				for line in in_reader:
				
					# skip track names and headers
					if len(line) == 1:
						track_name = get_track_name(line[0])
						continue
					
					if line[0] == '#chrom':
						continue
			
					# ucsc browser has a problem if end is less than 1
					if line[2] == '0':
						line[2] = '1'
				
					# add sample name at end of read IDs
					line[3] = line[3] + f"___{track_name}"
					out_writer.writerow(line)

		# sort integration sites so that we can collect the ones that are the same
		# between samples
		print("sorting integration sites")
		sp.run(['sort', '-k1,1', '-k2,2n', '-o', tmpfile, tmpfile])
		
		
		print("combining identical integration sites from different samples")
		with open(outfile, 'w') as outhandle, open(tmpfile, 'r') as tmphandle:
		
			# write header
			outhandle.write('track name=int_sites description="integration sites"\n')
			outhandle.write('#chrom\tChromStart\tChromEnd\tname\n')
			
			# create csv writer and reader objects
			in_reader = csv.reader(tmphandle, delimiter='\t')
			out_writer = csv.writer(outhandle, delimiter='\t')
			
			try:
				last = next(in_reader)
				reads = last[3]
			except StopIteration:
				print("didn't find any lines in tmpfile {tmpfile}!")
			
			combined = 0
			
			for line in in_reader:
			
				# if this site is different to the last one we saw
				# write it
				if line[:3] != last[:3]:
					last = line

					last[3] = reads
					out_writer.writerow(last)
					reads = line[3]
					
				# if it's the same
				else:
					reads = reads + f"+++{line[3]}"
					combined += 1
					
			# write last line
			last[3] = reads
			out_writer.writerow(last)	
			print(f"combined {combined} sites")
			print(f"saved output to {outfile}")	
					
					
	sp.run(['rm', tmpfile])

 
def get_track_name(line):
	track_name = re.match("track name=(.+) description=", line)
	try:
		track_name = track_name.group(1)
	except IndexError:
		print(f"couldn't find track name in line {line}")
		track_name = ''	
		
	return track_name
	
	
if __name__ == "__main__":
	main()

