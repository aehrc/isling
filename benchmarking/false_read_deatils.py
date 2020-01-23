### finds the seq, CIGAR, XA and SA of reads ##

## pipeline appears to produce more false positives and false negatives for reads containing integrations with overlaps 
## to determine whether issues are a result of error in the intial insert_virus.py or read_integrations.py scripts can can further 
## analyse these false reads using the corresponding seq, CIGAR, XA (secondary alignment) and SA (supplementary alignment) 

##import libraries 
import pandas as pd
from simplesam import Reader
import argparse
import sys


##main
def main(argv): 

	print("STARTING...")
	parser = argparse.ArgumentParser(description='obtain seq, CIGAR, XA and SA for a list of read IDs') 
	parser.add_argument('--ID', help = 'list of read IDs missed by the pipeline', required = True)
	parser.add_argument('--bam', help = 'host aligned bam file', required = True)
	args = parser.parse_args() 


	#get the list of read IDs 
	false_ID = open(args.ID, 'r') 
	false_ID = false_ID.read().splitlines()


	#get details on the reads 
	false_details = getDetails(args.bam, false_ID) 

	#save dataframe of false details as a dataframe 
	false_details.to_csv('read_details.csv', sep = '\t') 
	print('Details of parsed reads saved to read_details.csv') 

def getDetails(bam, false_IDs): 
	"""Extract the required details from the bam file on each of the missed integrations""" 
		
	#create reader to read through bam file
	in_file = open(bam, 'r') 
	in_bam = Reader(in_file)

	#get the number of reads in the bam file
	num_bam = len(in_bam) 

	#lists for the information which we ned 
	read_id = []
	read_seq = [] #sequence of the false reads
	read_pos = [] #leftmost position of read
	read_cig = [] #cigar for the false reads
	read_xa = [] #list of secondary alignment
	read_sa = [] #list of primary alignment 


	#loop through the reads in the bam file 
	for i in range(num_bam): 
		x = next(in_bam) 
		if x.qname in false_IDs: 
			read_id.append(x.qname) 
			read_seq.append(x.seq)
			read_pos.append(x.pos) 
			read_cig.append(x.cigar)
			xa, sa = processTags(x.tags)  
	
	#create dataFrame of information
	false_details = pd.DataFrame({"ID": read_id,"Pos": read_pos, "Seq":read_seq, "Cigar": read_cig, "XA": xa, "SA": sa}) 
	
	return false_details

def processTags(tags): 
	"""Breaks up the tags for a read to obtain XA and SA""" 
	
	#get tags 
	xa = tags.get('XA')
	sa = tags.get('SA')
 
	return xa, sa



if __name__ == "__main__":
	main(sys.argv[1:]) 




