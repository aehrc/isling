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

	num_reads = numReads(read_file) 
	int_coord = intCoords(int_file)
	sam_file = args.sam
	read_id, left_read, right_read = processReads(sam_file,num_reads) 
	left_status, right_status, int_type, left_viral, right_viral = findOverlaps(left_read, right_read, int_coord)
	
	#save the file 
	results = pd.DataFrame({"read_id":read_id,"left_read_virus":left_status,"right_read_virus":right_status,"integration_type":int_type,"viral_bases_left":left_viral, "viral_bases_right":right_viral})
	with open(args.save, 'w') as handle: 
		results.to_csv(handle,sep='\t') 
		
	print("COMPLETE")
	print("Saved to "+str(args.save))
		
def numReads(sam_file): 
	"""Finds number of reads in the .sam file"""
	df = pd.read_csv(sam_file,header=None) 
	num_reads = len(df[~df[0].str.contains('@')])
	num_inserts = int(num_reads/2)
	return num_inserts 
	
def findOverlaps(left_read,right_read, int_coord): #try remaking this entire function 
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
			check_left = checkOverlap(left_read[i],int_coord[j])
			if check_left == True:
				overlap 
			
			
			check_right = checkOverlap(right_read[i], int_coord[j]) 
			if check_left == True or check_right == True:
				status = check_overlap
				#important if a read has more than one read (rare but possible) 
				overlap = overlap + overlapLength(read_coord[i],int_coord[j])
				
		overlap_status.append(status)
		overlap_len.append(overlap)
		
	return overlap_status, overlap_len
	
def findOverlaps(left_read, right_read, int_coord): 
	"""finds whether a read contains viral DNA (True/False) and the amount of viral DNA in the read (bp)"""
	#list of overlap statuses 
	left_int = []
	right_int = []
	
	#list of the amount of viral DNA in each read 
	left_viral = []
	right_viral = []
	
	#list of read types
	int_type = []
	
	for i in range(len(left_read)): 
		check_left = False
		check_right = False 
		l_coord = (-1,-1)
		r_coord = (-1,-1)
		l_overlap = 0 
		r_overlap = 0 
		for j in range(len(int_coord)):
		
			#check left read for viral DNA 
			c_left = checkOverlap(left_read[i],int_coord[j])
			if c_left == True: 
				check_left = True
				l_coord = j #the integetation which causes overlap on the left
				 
				#if overlaps find the length of said overlap 
				l_overlap = overlapLength(left_read[i], int_coord[j]) 
				
			#check right read for viral DNA
			c_right = checkOverlap(right_read[i], int_coord[j])
			if c_right == True: 
				check_right = True
				r_coord = j #the integration which causes overlap on the right 
				
				#if overlaps find the length of said overlap 
				r_overlap = overlapLength(right_read[i],int_coord[j])
				

		left_int.append(check_left)
		right_int.append(check_right) 	
		read_type = readType(check_left, l_coord, check_right, r_coord)
		int_type.append(read_type)
		left_viral.append(l_overlap)
		right_viral.append(r_overlap)
		
	return left_int, right_int, int_type, left_viral, right_viral 
		   

			
def readType(check_left, l_coord, check_right, r_coord):
	"""Finds the type of overlap ie - chimeric, split end end or all viral""" 
	type = ""
	
	#check if read is chimeric 
	if check_left == True and check_right == False or check_left == False and check_right == True: 
		type = "chimeric"
		
	elif check_left == True and check_right == True: 
		#if viral DNA spans the read
		if l_coord == r_coord: 
			type = "all viral DNA"
		#consider two different viral integrations at both ends 
		#this is unlikely but we consder for completeness 
		else: 
			type = "split ends"   
	
	elif check_left == False and check_right == False: 
		type = "no viral DNA"
	
	return type
		
def intCoords(int_file): 
	"""Finds the location of viral DNA in the host sequence""" 
	int_coord = []

	for i in range(len(int_file)):
		c1 = int_file["Start point"][i]
		c2 = int_file["Stop point"][i]
		int_coord.append((c1,c2))
	return int_coord
	
def processReads(sam_file,num_inserts):
	"""Creates a list of the read IDs and their alignment coordinates""" 

	read_file = open(sam_file,'r')
	in_sam = Reader(read_file) 
	#list of insert IDs
	read_id = []

	#list of insert coordinates in the original fasta sequence
	left_read = []
	right_read = []
	for i in range(0,num_inserts):
		if i%500000==0 and i!=0:
			print(str(i)+" READS PROCESSED") 
		x = next(in_sam)
		y = next(in_sam)
		
		#save the ID of the read 
		read_id.append(x.qname)
		
		#save the coordinates of the read 
		left_read.append((x.pos-1,x.pos+x.tlen-1))
		right_read.append((y.pos-1,y.pos+y.tlen-1))
		
	return read_id, left_read, right_read
	

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
	
	#find the lower point of A2 and B2 
	upper = min(B2, A2) 
	
	#find the upper point of A1 and B1 
	lower = max(B1, A1) 
	
	overlap = upper-lower 
	if overlap<0: 
		print("coordinate a (read): "+str(coordA))
		print("coordinate b (integration): "+str(coordB))# for debugging remove later 
		raise OSError("Attempting to find length of overlap between regions which do not overlap")

	return overlap
	
	
if __name__ == "__main__":
	main(sys.argv[1:]) 	
