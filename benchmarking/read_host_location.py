#converts location in host to location in host genome 
#can be added into one of the other scripts once finished  

#import libraries 
import pandas as pd
import numpy as np 
import argparse
import sys 
import csv


##main
def main(argv): 

	#get arguments
	parser = argparse.ArgumentParser(description='the location of reads in the host alignment')
	parser.add_argument('--reads', help = 'csv containing information the integrations in reads', required = True)
	parser.add_argument('--host_ints', help = 'file containing locations of host integrations', required = True) 
	args = parser.parse_args()


	#make dataframe from host ints file
	reads = pd.read_csv(args.reads, header = 0 , sep='\t') 

	#read in the host ints file 
	host_ints = pd.read_csv(args.host_ints, header = 0 , sep='\t') 	

	#make a dictionary of the length of each integration 
	int_lens = [host_ints['Stop point'][i]-host_ints['Start point'][i] +1 for i in range(len(host_ints))]
	int_dict = pd.Series(int_lens,index = host_ints['Start point'])

	#adjust right reads
	adj_Rreads = []


	for i in range(len(reads)): 

		#process the coordinates of all reads except for viral reads as they can not be mapped to the host 
		if reads["right_read"][i] != 'v': 

			#reads in the reverse orientation must be handled separately 
			if reads["left_read"][i] != 'v' ad reads["right_read"][i] != 'hv': 
				left_coor = reads['left_read_coor'][i].split(", ")
				#obtain the intial start and stop positions before adjustment 
				intStart = int(left_coor[1].replace('(',''))
				intStop = int(right_coor[1].replace('(',''))	

			else: 
				right_coor = reads['right_read_coor'][i].split(", ")

				#obtain the intial start and stop positions before adjustment 
				intStart = int(right_coor[0].replace('(','')) 
				intStop = int(right_coor[1].replace(')','')) 
		
			#loop through integrations to get location of reads in host 
			start = intStart 
			stop = intStop 

			for ints in host_ints['Start point']:
				if ints < intStart: 
					start = start - int_dict.get(ints)
					stop = stop - int_dict.get(ints)

			#adjust the start and stop positions as the SAM/BAM file format uses a 1-based coordinate system
			start = start+1 
			stop = stop+1 
			
		
		#give some placeholder values to viral rads 
		else: 
			start, stop = (-1,-1) 

	
		#store the adjusted coordinates 		
		adj_Rreads.append((start,stop)) 
	

	#adjust left reads 
	adj_Lreads = []

	for i in range(len(reads)): 

		#process the coordinates of all reads except for viral reads as they can not be mapped to the host
		if reads["left_read"][i] != 'v': 
			left_coor = reads['left_read_coor'][i].split(", ")

			#obtain the intial start and stop positions before adjustment 
			intStart = int(left_coor[0].replace('(','')) 
			intStop = int(left_coor[1].replace(')','')) 

			#loop through integrations to get the location of reads in the host 
			start = intStart 
			stop = intStop 

			for ints in host_ints['Start point']:
				if ints < intStart: 
					start = start - int_dict.get(ints)
					stop = stop - int_dict.get(ints)

			#adjust the start and stop positions as the SAM/BAM file format uses a 1-based coordinate system 
			start = start +1 
			stop = stop +1 

		else: 
			start, stop = (-1,-1) #placeholder for viral reads

		#store the adjusted coordinates 		
		adj_Lreads.append((start,stop)) 

	
	#save the adjusted coordinates as a dataframe 
	reads["host_left"] = pd.Series(adj_Lreads,index = reads.index)
	reads["host_right"] = pd.Series(adj_Rreads,index = reads.index) 

	#replace placeholder -1's for viral reads with nan
	#reads =  reads.assign(host_left = reads.host_left.where(reads.host_left.ge(0))) 
	#reads = reads.assign(host_right = reads.host_right.where(reads.host_right.ge(0)))

	#save to csv file 
	reads.to_csv('host_location_reads.csv', sep = '\t') 

	print("COMPLETE") 

if __name__ == "__main__":
	main(sys.argv[1:]) 
