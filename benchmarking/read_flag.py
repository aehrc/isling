#### assesses reads for the corresponding flag 

### import libraries
from simplesam import Reader
import pandas as pd
import argparse 
import sys

### def main
def main(argv): 
	print("Starting...") 
	
	#get arguments
	parser = argparse.ArgumentParser(description='determine the flags of the missed reads') 
	parser.add_argument('--bam', help = 'bam file of host aligned reads', required = True)
	parser.add_argument('--missed_reads', help = 'csv file containing the IDs of the reads which the pipeline failed to identify', required = True)
	parser.add_argument('--all_reads', help = 'csv file containing all reads', required = False)

	args = parser.parse_args()

	#bam file 
	bam_file = args.bam 

	#procsess reads 
	processReads(bam_file) 

	print("DONE")


def processReads(bam_file):
	"""make a list of lists of the flags for each read ID""" 

	#reader for the bam file 
	read_file = open(bam_file,'r')
	in_bam = Reader(read_file) 	

	#number of entries in bam file 
	entries = len(in_bam)

	#list of IDs 
	read_ID = []

	#list of flags 
	read_flag = []

	#list of pos 
	read_pos = []  

	# loop through the bam file entries 
	for i in range(entries): 
		
		#report the number of reads processed 
		if i != 0 and i%100000 == 0: 
			print(str(i) + " reads processed") 

		x = next(in_bam) 
		read_ID.append(x.qname) 
		read_flag.append(x.flag)
		read_pos.append(x.pos -1 )  #-1 to match output from read_integrations.py 
 
	#make dataframe of IDs and flags 	
	all_bam = pd.DataFrame({"ID":read_ID, "Flag":read_flag,"Pos":read_pos}) 
	
	#remove the flags corresponding to 2048 - "supplementary alignment" and 2064 - "supplementary alignment -  read reverse strand" 
	drop_flag = [] #flags to be dropped 

	for i in range(len(all_bam)): 
		if all_bam["Flag"][i] == 2048 or all_bam["Flag"][i] == 2064: 
			drop_flag.append(i) 
 
	all_bam = all_bam.drop(drop_flag)

 	#create a list of unique IDs 
	unique_ID = all_bam["ID"].drop_duplicates()
	unique_ID = unique_ID.reset_index(drop=True) 
	all_bam.set_index(["ID"], inplace = True, drop = True)

	return all_bam 
	

if __name__ == "__main__":
	main(sys.argv[1:]) 
