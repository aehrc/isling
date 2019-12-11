### combines data from simiulated reads and integrated host sequence to output which reads contain viral DNA ###
## uses .sam file outputted from reads simulated using ART https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3278762/
## uses file output from insert_virus.py which lists the location of viral DNA in the human .fasta


#import libraries 
from simplesam import Reader 
import pandas as pd
import argparse
import sys


### main
def main(argv):

	#get arguments 
	parser = argparse.ArgumentParser(description='determine reads which contain viral DNA')
	parser.add_argument('--sam', help = 'sam file describing simulated reads', required = True)
	parser.add_argument('--host_ints', help = 'csv file describing the location of viral DNA in the human host', required = True)
	parser.add_argument('--save', help = 'csv file to save the output', required = True)
	args = parser.parse_args() 
	
	#read in sam file
	read_file = open(args.sam,'r')
	in_sam = Reader(read_file) 
	
	#read in host integration locations
	int_file = pd.read_csv(args.host_ints,header=0,sep='\t')
	
	print("DOING TASKS")
	num_reads = numReads(read_file) 
	print('TASK ONE')
	int_coord = intCoords(int_file)
	print('TASK TWO')
	read_id, read_coord = processReads(in_sam,num_reads) 
	print('TASK THREE')
	overlap_status, overlap_len = findOverlaps(read_coord,int_coord)
	
	#save the file 
	results = pd.DataFrame({"read_id":read_id,"contains_virus":overlap_status,"virus_length":overlap_len})
	#handle = open(args.locs,"w+")
	#handle.write("HI")
	print(results[0:3])
	results.to_csv("read_integrations.csv",sep='\t')
	#print("THING SAVING TO "+args.locs)
	
	
def numReads(sam_file): 
	"""Finds number of reads in the .sam file"""
	df = pd.read_csv(sam_file,header=None) 
	num_reads = len(df[~df[0].str.contains('@')])
	num_inserts = int(num_reads/2)
	return num_inserts 
	
def findOverlaps(read_coord,int_coord): 
	"""Creates a list of whether a read overlaps (True/False) and a list of the length of the corresponding overlap
	takes lists of the coordinates of the reads and the locations of the viral DNA in the host"""
	#list of overlap statuses 
	overlap_status = []

	#list of overlap lengths 
	overlap_len = []
	
	for i in range(len(read_coord)):
		overlap = 0
		status = False
		for j in range(len(int_coord)):
			status = checkOverlap(read_coord[i],int_coord[j])
			if status == True:
				overlap = overlapLength(read_coord[i],int_coord[j])
				break
		overlap_status.append(status)
		overlap_len.append(overlap)
		
	return overlap_status, overlap_len	

def intCoords(int_file): 
	"""Finds the location of viral DNA in the host sequence""" 
	int_coord = []

	for i in range(len(int_file)):
		c1 = int_file["Start point"][i]
		c2 = int_file["Stop point"][i]
		int_coord.append((c1,c2))
	return int_coord
	
def processReads(in_sam,num_inserts):
	"""Creates a list of the read IDs and their alignment coordinates""" 

	#list of insert IDs
	read_id = []

	#list of insert coordinates in the original fasta sequence
	read_coord = []
	print("NUM INSERTS "+str(num_inserts))
	for i in range(0,5): 
		x = next(in_sam)
		y = next(in_sam)
		read_id.append(x.qname)
		read_coord.append((x.pos-1,x.pos+x.tlen-1))
		
	return read_id, read_coord
	

def checkOverlap(coordA,coordB):
	"""Tells us whether coordA overlaps with coordB"""
	status = False 
	A1, A2 = coordA
	B1, B2 = coordB 
	if A2>B1 and A1<B1:
		status = True 
	return status 

def overlapLength(coordA,coordB):
	"""Tells us the length of the overlap betweeen coordA and coordB"""
	A1, A2 = coordA
	B1, B2 = coordB
	overlap = 0
	#isolate where overlaps occur
	if A2>B1 and A1<B1:
        	#if overlap is at the left end of B
		if A1<B1: 
			overlap = A2-B1
		#if overlap is at the right end of B
		if A2>B2: 
			overlap = B2-A1
		#if overlap is within B 
		else: 
			overlap = A2-A1
	return overlap
	
if __name__ == "__main__":
	main(sys.argv[1:]) 	
