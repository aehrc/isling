### combines data from simiulated reads and integrated host sequence to output which reads contain viral DNA ###

## uses .sam file outputted from reads simulated using ART https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3278762/
## uses file output from insert_virus.py which lists the location of viral DNA in the human .fasta
## when using 'coordinates' we consider (start point, end point) of either the viral integration in the host genome or the reads when aligned 


#import libraries 
from simplesam import Reader 
import pandas as pd
import argparse
import sys
import csv


### main
def main(argv):

	print("STARTING...", flush=True)

	#get arguments 
	parser = argparse.ArgumentParser(description='determine reads which contain viral DNA')
	parser.add_argument('--sam', help = 'sam file describing simulated reads', required = True)
	parser.add_argument('--host_ints', help = 'csv file describing the location of viral DNA in the human host', required = True)
	parser.add_argument('--save', help = 'csv file to save the output', required = True)
	parser.add_argument('--viral_only', help = 'csv file to save only information for viral reads', required = True)
	parser.add_argument('--filtered', help = 'list of IDs which are not mapped', required = False, default = False) 
	args = parser.parse_args()
	
	#read in sam file
	read_file = open(args.sam,'r')
	in_sam = Reader(read_file)
	
	#read in host integration locations
	#read in as chunks as the file may be very large  
	chunks = pd.read_csv(args.host_ints,header=0,sep='\t', chunksize = 50000) 
	int_file = pd.concat(chunks) 

	#read in the list of filtered read IDs 
	if args.filtered != False: 
		id_file = open(args.filtered, 'r')
		filtered_id = id_file.read().splitlines()
	else: 
		filtered_id = []


	#find the number of fragments 
	num_reads = numReads(read_file)

	#get the coordinates of each integration
	int_coord = intCoords(int_file)

	#get the hPos of each integration - we use this later as an identifer
	int_hPos = intHpos(int_file)

	#get the types of junctions for each integration
	int_leftj, int_rightj = intJunction(int_file) 
	
	#read in sam file to process reads 
	sam_file = args.sam

	#process reads to obtain the ID of each fragment and the coordinates of the corresponding reads
	fragment_id, first_read, second_read = processReads(sam_file,num_reads, filtered_id)

	#assess for the type of read and the amount of viral DNA (bp) in each read
	first_type, second_type, first_len, second_len, first_loc, second_loc, read_hPos, first_junc, second_junc  = analyseRead(first_read,second_read, int_coord, int_hPos, int_leftj, int_rightj)
	

	#save the entire file
	results = pd.DataFrame({"fragment_id":fragment_id,"left_read":first_type,"right_read":second_type,"left_read_amount":first_len,"right_read_amount":second_len,
"left_read_coor":first_read,"right_read_coor":second_read, "left_read_Vcoor": first_loc, "right_read_Vcoor": second_loc, "hPos": read_hPos, "left_junc":first_junc, "right_junc":second_junc})
	with open(args.save, 'w') as handle: 
		results.to_csv(handle,sep='\t')

	#create a file which saves information only for reads which contain viral DNA
	viral_only = getViralReads(results)
	with open(args.viral_only, 'w') as handle: 
		viral_only.to_csv(handle, sep='\t')
		
	print("COMPLETE")
	print("Information on ALL reads saved to: "+str(args.save))
	print("Information on VIRAL reads saved to: "+str(args.viral_only))
	
		
def numReads(sam_file): 
	"""Finds number of reads in the .sam file"""

	df = pd.read_csv(sam_file,header=None) 
	num_reads = len(df[~df[0].str.contains('@')])
	num_inserts = int(num_reads/2)
	return num_inserts 
	
		
def analyseRead(first_read,second_read, int_coord, int_hPos, int_leftj, int_rightj):  
	"""Creates a list of whether a read overlaps (True/False) and a list of the length of the corresponding overlap
	takes lists of the coordinates of the reads and the locations of the viral DNA in the host"""
	# this function is very intensive - requires optimisiing 	

	print("\nANALYSING READS", flush = True) 

	#lists the types of the reads  
	first_type = []
	second_type = []
	
	#lists the amount of viral DNA in each read
	first_len = []
	second_len = []
	
	#list the location of an integration in the read 
	first_loc = []
	second_loc = []
	
	#list of the integrations in each read - denoted as hPos 
	read_hPos = []

	#list the junction types in each read
	first_junc = []
	second_junc = [] 
	
	#loop to iterate through the reads. len(first_read) used though could have used len(second_read) 
	for i in range(len(first_read)):
		#count number of bases which overlap 
		overlap_len1 = []
		overlap_len2 = []
		
		#keep record of which overlaps occur
		overlap_type1 = []
		overlap_type2 = []

		#record the location of viral DNA in reads
		readloc_1 = []
		readloc_2 = []

		#list the integrations in the read 
		int_list = []

		#list the junction types 
		junction_1 = ""
		junction_2 = ""
		
		for j in range(len(int_coord)):
			
			#compare the ith left read with the jth integration 
			c_type = checkOverlap(first_read[i],int_coord[j])			
			if c_type != "": 
				#store information on the type of integration 
				overlap_type1.append(c_type) 
				
				#save information on the amount of viral DNA in the read from the integration 
				overlap, start_1, stop_1 = overlapLength(first_read[i],int_coord[j])
				overlap_len1.append(overlap)
				readloc_1.append((start_1, stop_1)) 

				#Store the hPos of the integration in the read 
				int_list.append(int_hPos[j])  
			
			#compare the ith right read with the jth integration 
			c_type = checkOverlap(second_read[i], int_coord[j]) 
			if c_type != "": 
				#store information on the type of integration 
				overlap_type2.append(c_type) 
				
				#save information on the amount of viral DNA in the read from the integration 
				overlap, start_2, stop_2 = overlapLength(second_read[i], int_coord[j]) 
				overlap_len2.append(overlap)
				readloc_2.append((start_2,stop_2))

				#Store the hPos of the integration in the read 
				int_list.append(int_hPos[j]) 
		
		#find the types of the read 
		type1 = readType(overlap_type1)
		type2 = readType(overlap_type2) 
		
		#save these read types 
		first_type.append(type1) 
		second_type.append(type2)
		
		#find the amount of viral DNA in each read 
		len1 = viralQuantity(overlap_type1, overlap_len1)
		len2 = viralQuantity(overlap_type2, overlap_len2) 
		
		#save the amount of viral DNA in each read 
		first_len.append(len1)
		second_len.append(len2)

		#save the start and stop positions of viral DNA in each read 
		first_loc.append(readloc_1)
		second_loc.append(readloc_2) 

		#save the integrations (hPos) in each read 
		read_hPos.append(set(int_list)) #exists as a list of lists

		#save the junction types of each read
		junction_1 = readJunction(j, overlap_type1 ,int_leftj, int_rightj)
		junction_2 = readJunction(j, overlap_type2 ,int_leftj, int_rightj) 
		first_junc.append(junction_1)
		second_junc.append(junction_2)

		#report how many reads have been analysed 
		if i % 250000 == 0: 
			print("{:.2f}".format((i*100/len(first_read)))+"% of reads analysed...",flush = True) 	
		
	return first_type, second_type, first_len, second_len, first_loc, second_loc, read_hPos, first_junc, second_junc 

def readType(overlap_type): 
	
	#intialise string
	read_type = ""
	
	#handle reads which are all viral DNA 
	if "all" in overlap_type: 
		read_type = "v" 

	#handle short reads - where there is viral DNA in the middle of a read 	
	elif "short" in overlap_type: 
		read_type = "sh" 
	
	#handle chimeric reads
	elif "left" in overlap_type and "right" not in overlap_type: 
		read_type = "vh"

	#handle chimeric reads 
	elif "right" in overlap_type and "left" not in overlap_type: 
		read_type = "hv" 
	
	#handle reads without viral DNA
	else: 
		read_type = "h" 
	
	return read_type
	
		
def viralQuantity(overlap_type, overlap_len): 
	"""function which finds the amount of viral DNA (bp) in a read""" 
	
	#intalise amount of viral DNA 
	viral_q = 0
	
	#handle reads which span an integration
	if "all" in overlap_type: 
		viral_q = max(overlap_len) 		
 		
	#handle chimeric, split end or short reads
	#it is possible (though extremely unlikely) a read could be short and have split ends 
	#don't consider this in the other function as short read is more important than split ends but consider for when we calculate the number of viral bp in  a read 
	else: 
		viral_q = sum(overlap_len) 
		
	return viral_q		   
		
def intCoords(int_file): 
	"""Finds the location of viral DNA in the host sequence"""
	#get the integration coordinates  
	int_coord = []

	for i in range(len(int_file)):
		c1 = int_file["Start point"][i] 
		c2 = int_file["Stop point"][i]
		int_coord.append((c1,c2))

	return int_coord

def intHpos(int_file): 
	"""finds the hPos of the host integration allowing subsequent identification""" 

	#list of hPos values for each integration in int_coord[]
	int_hPos = []

	for i in range(len(int_file)): 
		hPos = int_file['hPos'][i]
		int_hPos.append(hPos) 

	return int_hPos

def intJunction(int_file):
	"""finds the left and right junction types for each integration""" 
	
	#list of the junctions 
	int_leftj = []
	int_rightj = []

	for i in range(len(int_file)):
		#get junctions  
		leftj = int_file['leftJunction'][i]
		rightj = int_file['rightJunction'][i]
		#save the junctions
		int_leftj.append(leftj) 
		int_rightj.append(rightj)

	return int_leftj, int_rightj

def readJunction(index, overlap_type, int_leftj, int_rightj): 
	"""finds the type of junction in a read""" 

	#if entire stretch is viral DNA there is no junction
	if "all" in overlap_type: 
		junction = "none" 

	#if a read is a short read it has two junctions which is complicated 
	if "short" in overlap_type: 
		junction = "other" 
	
	#if viral DNA is on the right we care about the left junction
	elif "right" in overlap_type and "left" not in overlap_type: 
		junction = int_leftj[index]

	#if viral DNA is on the left we care about the right junction 
	elif "left" in overlap_type and "right" not in overlap_type: 
		junction =  int_rightj[index]

	#it is complicated to handle split reads as there is more than one junction  
	elif "left" in overlap_type and "right" in overlap_type: 
		junction = "other"

	#otherwise we assume that we have host DNA only in the fragment 
	else: 
		junction = ""

	return junction

	
def processReads(sam_file,num_inserts, filtered_id):
	"""Creates a list of the read IDs and their alignment coordinates""" 


	read_file = open(sam_file,'r')
	in_sam = Reader(read_file) 
	#list of insert IDs
	fragment_id = []

	#list of insert coordinates in the original fasta sequence
	first_read = []
	second_read = []


	print("PROCESSING READS", flush = True) 
	for i in range(0,num_inserts):
		if i%500000==0 and i!=0:
			print("{:.2f}".format((i*100/num_inserts))+"% of reads processed...") 
		x = next(in_sam)
		y = next(in_sam)

		#get ID of the read being processed
		this_ID = str(x.qname) 
		
		#only consider the read is mapped 
		if len(filtered_id) > 0: 
			if this_ID in filtered_id: 
				#save the ID of the read 
				fragment_id.append(x.qname)
		
				#save the coordinates of the read 
				#subract 1 as the SAM file position starts at 1 but we use index 0 
				first_read.append((x.pos-1,x.pos+len(x.seq)-1))
				second_read.append((y.pos-1,y.pos+len(y.seq)-1))
		else: 
			#save the ID of the read 
			fragment_id.append(x.qname)
		
			#save the coordinates of the read
			#subtract 1 as the SAM file position starts at 1 but we use index 0  
			first_read.append((x.pos-1,x.pos-1+(len(x.seq)-1)))
			second_read.append((y.pos-1,y.pos-1+(len(y.seq)-1)))
			
	return fragment_id, first_read, second_read
	

def checkOverlap(coordA, coordB): 
	""" Function which tells us whether coordA (start of read, end of read) overlaps coordB (start of integration, end of integrations). Returns string which says which type of integrations occured""" 
	
	#get coordinates 
	A1, A2 = coordA #start and stop position of read 
	B1, B2 = coordB #start and stop position of integration
	

	if B1<=A1 and B2>=A2: 
		overlap_type = "all"    #virus covers the entire read 
	elif A1 < B1 and A2 > B2: 	 
		overlap_type = "short"  #virus is in the centre of the read 
	elif A1 < B1 and B1 <= A2:
		overlap_type = "right" 	#virus is on the right end of the read 
	elif B2 < A2 and B2 >= A1:
		overlap_type = "left"   #virus on the left end of the read 
	else:
		overlap_type = ""
		
	return overlap_type 


def getViralReads(results): 
	"""Function to create dataframe consisting of only viral reads. While this is more coded than required, this is the most efficient way to create a large dataframe with the required information""" 
	
	#create a list of the indexes we want to drop (more code but more efficient) 
	idx = []
	
	#reindex column 
	results = results.reset_index(drop=True)

	for i in range(len(results)): 
		#no integration if left and right reads are host
		if results['left_read'][i] == 'h' and results['right_read'][i] == 'h': 
			idx.append(i)
		#no integration if left and right reads are viral
		if results['left_read'][i] == 'v' and results['right_read'][i] == 'v':
			idx.append(i)

	viral_reads = results.drop(results.index[idx]) 
	print("Number of viral reads: "+str(len(viral_reads)))

	return viral_reads

		 		
def overlapLength(coordA,coordB): 
	"""Tells us the length of the overlap betweeen coordA and coordB"""	

	A1, A2 = coordA #start and stop position of read
	B1, B2 = coordB #start and stop position of integration 

	#find the amount of DNA shared by the read and the integration
	#plus 1 to find the number of viruses (ie if the first base of an integration was at 4 and the last base at 6 the number of viral bases is 6 - 4 +1 = 3 as bases 4,5 and 6 are viral.  
	overlap = min(B2, A2)-max(B1, A1)+1

	if overlap<0: 
		raise OSError("Attempting to find length of overlap between regions which do not overlap")

	#if there is an integration find the position of the integration in the read 
	vStart = ""
	vStop = ""
	if overlap > 0: 
		vStart, vStop = vStartStop(coordA, coordB) 
	
	return overlap, vStart, vStop

def vStartStop(coordA, coordB): 
	"""Finds the start and stop positions of viral DNA in a read"""
	
	A1, A2 = coordA #start and stop position of read
	B1, B2 = coordB #start and stop position of integration

	#get start position of integration in read 
	if A1 < B1:
		vStart = B1
	else: 
		vStart = B2

	#get stop position of integration in read 
	if A2 > B2: 
		vStop = B2 
	else: 
		vStop = A2 

	#adjust start and stop position so that it is given in relation to the read
	vStart = vStart - A1
	vStop = vStop - A1
	
	return vStart, vStop 

if __name__ == "__main__":
	main(sys.argv[1:]) 	
