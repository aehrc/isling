### combines data from simiulated reads and integrated host sequence to output which reads contain viral DNA ###
## uses .sam file outputted from reads simulated using ART https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3278762/
## uses file output from insert_virus.py which lists the location of viral DNA in the human .fasta
## when using 'coordinates' we consider (start point, end point) of either the viral integration in the host genome or the reads when aligned 


#import libraries 
from simplesam import Reader 
import pandas as pd
import argparse
import sys
import time #remove after debugging - used to optimise code TODO


### main
def main(argv):

	print("STARTING", flush=True)

	#get arguments 
	start_time = time.time()
	parser = argparse.ArgumentParser(description='determine reads which contain viral DNA')
	parser.add_argument('--sam', help = 'sam file describing simulated reads', required = True)
	parser.add_argument('--host_ints', help = 'csv file describing the location of viral DNA in the human host', required = True)
	parser.add_argument('--save', help = 'csv file to save the output', required = True)
	parser.add_argument('--viral_only', help = 'csv file to save only information for viral reads', required = True) 
	args = parser.parse_args()
	now = time.time()
	print("Time taken to parse arguments "+format(now-start_time),flush=True) 
	
	#read in sam file
	read_file = open(args.sam,'r')
	in_sam = Reader(read_file)
	now = time.time() 
	print("Time taken to read in sam file "+format(now-start_time),flush=True)
	
	#read in host integration locations
	start_time = time.time()
	#read in as chunks as the file may be very large - TODO optimise this 
	chunks = pd.read_csv(args.host_ints,header=0,sep='\t', chunksize = 50000)
	int_file = pd.concat(chunks)  
	now = time.time() 
	print("Time taken to read in host integrations "+format(now-start_time),flush=True)


	#find the number of fragments 
	start_time = time.time()
	num_reads = numReads(read_file)
	now = time.time() 
	print("Time taken to find the number of fragments "+format(now-start_time),flush=True)

	#get the coordinates of each integration
	start_time = time.time() 
	int_coord = intCoords(int_file)
	now = time.time() 
	print("Time taken to find coordinates of each integration "+format(now-start_time),flush=True)

	#get the hPos of each integration - we use this later as an identifer
	start_time = time.time() 
	int_hPos = intHpos(int_file)
	now = time.time() 
	print("Time taken to find hPos of each integration "+format(now-start_time),flush=True) 

	#get the types of junctions for each integration
	start_time = time.time()
	int_leftj, int_rightj = intJunction(int_file)
	now = time.time() 
	print("Time taken to find the type of junction for each integration "+format(now-start_time),flush=True) 
	
	#read in sam file to process reads 
	start_time = time.time()
	sam_file = args.sam
	now = time.time() 
	print("Time taken to read in sam file "+format(now-start_time),flush=True) 

	#process reads to obtain the ID of each fragment and the coordinates of the corresponding reads
	start_time = time.time()
	fragment_id, first_read, second_read = processReads(sam_file,num_reads)
	now = time.time() 
	print("Time taken to process reads "+format(now-start_time),flush=True) 

	#assess for the type of read and the amount of viral DNA (bp) in each read
	start_time = time.time()
	first_type, second_type, first_len, second_len, read_hPos, first_junc, second_junc  = analyseRead(first_read,second_read, int_coord, int_hPos, int_leftj, int_rightj)
	now = time.time() 
	print("Time taken to analyse reads "+format(now-start_time),flush=True)
	
	#save the entire file
	start_time = time.time()
	results = pd.DataFrame({"fragment_id":fragment_id,"left_read":first_type,"right_read":second_type,"left_read_amount":first_len,"right_read_amount":second_len,
"left_read_coor":first_read,"right_read_coor":second_read, "hPos": read_hPos, "left_junc":first_junc, "right_junc":second_junc})
	with open(args.save, 'w') as handle: 
		results.to_csv(handle,sep='\t')
	now = time.time() 
	print("Time taken save entire file "+format(now-start_time),flush=True)
	
	#create a file which saves information only for reads which contain viral DNA
	start_time = time.time() 
	viral_only = getViralReads(results)
	with open(args.viral_only, 'w') as handle: 
		viral_only.to_csv(handle, sep='\t')
	now = time.time() 
	print("Time taken save viral only file "+format(now-start_time))
		
	print("COMPLETE")
	print("Information on ALL reads saved to: "+str(args.save))
	print("Information on VIRAL reads saved to: "+str(args.viral_only))
	print("Time taken to run entire program "+format(start_time)) 
	
		
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


	#lists of the types of the reads  
	first_type = []
	second_type = []
	

	#lists for the amount of viral DNA in each read
	first_len = []
	second_len = []
	
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
				overlap = overlapLength(first_read[i],int_coord[j])
				overlap_len1.append(overlap) 

				#Store the hPos of the integration in the read 
				int_list.append(int_hPos[j])

				#Store the junction type
				if junction_1 != "none":
					junction_1 = readJunction(j, c_type ,int_leftj, int_rightj)   
			
			#compare the ith right read with the jth integration 
			c_type = checkOverlap(second_read[i], int_coord[j]) 
			if c_type != "": 
				#store information on the type of integration 
				overlap_type2.append(c_type) 
				
				#save information on the amount of viral DNA in the read from the integration 
				overlap = overlapLength(second_read[i], int_coord[j]) 
				overlap_len2.append(overlap)

				#Store the hPos of the integration in the read 
				int_list.append(int_hPos[j])

				#Store the junction type
				if junction_2 != "none": 				
					junction_2 = readJunction(j, c_type ,int_leftj, int_rightj)  
		
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

		#save the integrations (hPos) in each read 
		read_hPos.append(int_list) #exists as a list of lists

		#save the junction types of each read
		first_junc.append(junction_1)
		second_junc.append(junction_2) 	
		
	return first_type, second_type, first_len, second_len, read_hPos, first_junc, second_junc 
	
	
def readType(overlap_type): 
	"""Identifies the type of a read: chimeric, split, viral or host.Uses one of the overlap_type lists (overlap_type1 or overlap_type2. Returns type as a string""" 
	
	#intialise string
	read_type = ""
	
	#handle reads which are all viral DNA 
	if "all" in overlap_type: 
		read_type = "viral" 

	#handle short reads - where there is viral DNA in the middle of a read 	
	elif "short" in overlap_type: 
		read_type = "short" 

	#handle chimeric reads 
	elif "left" in overlap_type and "right" not in overlap_type or "right" in overlap_type and "left" not in overlap_type: 
		read_type = "chimeric" 
		
	#handle split reads 
	#these are rare but we include them for completeness 
	elif "left" in overlap_type and "right" in overlap_type: 
		read_type = "split" 
		
	#hanlde reads without viral DNA 
	else: 
		read_type = "host" 
		
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

	junction = "" 
	#if viral DNA is on the right we care about the left junction
	if overlap_type == "right":
		junction = int_leftj[index]  	

	#if viral DNA is on the left we care about the right junction 
	if overlap_type == "left":
		junction = int_rightj[index] 

	else: 
		junction = "none" 

	return junction  
	
def processReads(sam_file,num_inserts):
	"""Creates a list of the read IDs and their alignment coordinates""" 

	read_file = open(sam_file,'r')
	in_sam = Reader(read_file) 
	#list of insert IDs
	fragment_id = []

	#list of insert coordinates in the original fasta sequence
	first_read = []
	second_read = []
	for i in range(0,num_inserts):
		if i%500000==0 and i!=0:
			print(str(i)+" READS PROCESSED") 
		x = next(in_sam)
		y = next(in_sam)
		
		#save the ID of the read 
		fragment_id.append(x.qname)
		
		#save the coordinates of the read 
		first_read.append((x.pos-1,x.pos+len(x.seq)))
		second_read.append((y.pos-1,y.pos+len(y.seq)))
		#print("COORD A: "+str(first_read))	
		#print("COORD B: "+str(second_read)) #coordinates apear in the wrong orientation 
	return fragment_id, first_read, second_read
	
 
	
def checkOverlap(coordA, coordB): 
	""" Function which tells us whether coordA (start of read, end of read) overlaps coordB (start of integration, end of integrations). Returns string which says which type of integrations occured""" 
	
	#intialise output string
	overlap_type = "" 
	
	#get coordinates 
	A1, A2 = coordA
	B1, B2 = coordB 
	

	if B1<=A1 and B2>=A2: 
		overlap_type = "all" 
	elif A1 < B1 and A2 > B2: 	 
		overlap_type = "short"  
	elif A1 < B1 and B1 <= A2:
		overlap_type = "right" 	#virus is on the right end of the read 
	elif B2 < A2 and B2 >= A1:
		overlap_type = "left"
		
	return overlap_type 
		 		

def overlapLength(coordA,coordB): 
	"""Tells us the length of the overlap betweeen coordA and coordB"""
	#debugging remove later 	

	A1, A2 = coordA
	B1, B2 = coordB
	
	#find the lower point of A2 and B2 
	upper = min(B2, A2) 
	
	#find the upper point of A1 and B1 
	lower = max(B1, A1) 
	
	overlap = upper-lower #subtract 1 to allow for difference resulting from coordinate system 
	if overlap<0: 
		print("coordinate a (read): "+str(coordA))
		print("coordinate b (integration): "+str(coordB))# for debugging remove later 
		raise OSError("Attempting to find length of overlap between regions which do not overlap")
	return overlap

def getViralReads(results): 
	"""Function to create dataframe consisting of only viral reads. While this is more coded than required, this is the most efficient way to create a large dataframe with the required information""" 
	
	#create a list of the indexes we want to drop (more code but more efficient) 
	idx = []
	
	#reindex column 
	results = results.reset_index(drop=True)

	for i in range(len(results)): 
		if results['left_read'][i] != 'host' or results['right_read'][i] != 'host': 
			idx.append(i)
	
	viral_reads = results.loc[idx] 
	print("Number of viral reads: "+str(len(viral_reads)))

	return viral_reads

	
if __name__ == "__main__":
	main(sys.argv[1:]) 	
