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
	fragment_id, first_read, second_read = processReads(sam_file,num_reads) 
	first_type, second_type, first_len, second_len = analyseRead(first_read,second_read, int_coord) 
	
	#save the file 
	results = pd.DataFrame({"fragment_id":fragment_id,"left_read":first_type,"right_read":second_type,"left_read_amount":first_len,"right_read_amount":second_len,})
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
	
	
def analyseRead(first_read,second_read, int_coord):  
	"""Creates a list of whether a read overlaps (True/False) and a list of the length of the corresponding overlap
	takes lists of the coordinates of the reads and the locations of the viral DNA in the host"""
	
	#lists of the types of the reads  
	first_type = []
	second_type = []
	

	#lists for the amount of viral DNA in each read
	first_len = []
	second_len = []
	
	
	#loop to iterate through the reads. len(first_read) used though could have used len(second_read) 
	for i in range(len(first_read)):
		#count number of bases which overlap 
		overlap_len1 = []
		overlap_len2 = []
		
		#keep record of which overlaps occur
		overlap_type1 = []
		overlap_type2 = []
		
		for j in range(len(int_coord)):
			
			#compare the ith left read with the jth integration 
			c_type = checkOverlap(first_read[i],int_coord[j])			
			if c_type != "": 
				#store information on the type of integration 
				overlap_type1.append(c_type) 
				
				#save information on the amount of viral DNA in the read from the integration 
				overlap = overlapLength(first_read[i],int_coord[j])
				overlap_len1.append(overlap) 
			
			#compare the ith right read with the jth integration 
			c_type = checkOverlap(second_read[i], int_coord[j]) 
			if c_type != "": 
				#store information on the type of integration 
				overlap_type2.append(c_type) 
				
				#save information on the amount of viral DNA in the read from the integration 
				overlap = overlapLength(second_read[i], int_coord[j]) 
				overlap_len2.append(overlap) 
		
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
		
	return first_type, second_type, first_len, second_len 
	
	
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
		first_read.append((x.pos-1,x.pos+x.tlen-1))
		second_read.append((y.pos-1,y.pos+y.tlen-1))
		
	return fragment_id, first_read, second_read
	
 
	
def checkOverlap(coordA, coordB): 
	""" Function which tells us whether coordA (start of read, end of read) overlaps coordB (start of integration, end of integrations). Returns string which says which type of integrations occured""" 
	
	#intialise output string
	overlap_type = "" 
	
	#get coordinates 
	A1, A2 = coordA
	B1, B2 = coordB 
	
	if A1 <= B1 and A2 >= B1:
		overlap_type = "left" 	
		if A2 >= B2: 
			overlap_type = "all"
	elif A1 <= B2 and A2>=B2:
		overlap_type = "right"
	elif A1 > B1 and A2 < B2: 
		overlap_type = "short"  
		
	return overlap_type 
		 
		

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
