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
	parser.add_argument('--locs', help = 'text file containing details of the integrations', required = True)
	parser.add_argument('--reads', help = 'csv containing information the integrations in reads', required = True)

	#make dataframe from locs file 
	locs = pd.read_csv(args.locs, header =0, sep='\t')

	#make dataframe from host ints file
	reads = pd.read_csv(args.reads, header = 0 , sep='\t') 

	
	#make a dictionary of the length of each integration 
	int_lens = [locs[vStop][i]-locs[vStart][i] for i in range(len(locs))]
	int_dict = pd.Series(int_lens,index = locs.vStart)

	
	#adjust right reads
	adj_Rread = []

	for i in range(len(reads)): 
		if reads["right_read"]!= 'v': 
			start, stop = read['right_read_coor'][i]
		else: 
			start, stop = (-1,-1) 
 
		for ints in locs['vStart']: 
			if ints < start: 
				start = start - int_dict.get(ints)
				stop = stop - int_dict.get(ints)
		adj_Rreads.append((start,stop)) 

	
	#adjust right reads 
	adj_Lread = [] 
	for i in range(len(reads)): 
		if reads["left_read"]!= 'v': 
			start, stop = read['left_read_coor'][i]
		else: 
			start, stop = (-1,-1) 
 
		for ints in locs['vStart']: 
			if ints < start: 
				start = start - int_dict.get(ints)
				stop = stop - int_dict.get(ints)
		adj_Lreads.append((start,stop))

	
	#save the adjusted coordinates as a dataframe 
	reads["host_left"] = pd.Series(adj_Lreads,index = reads.index)
	reads["host_right"] = pd.Series(adj_Rreads,index = reads.index) 

	#save to csv file 
	reads.to_csv('adjusted_reads.csv', sep = '\t') 


if __name__ == "__main__":
	main(sys.argv[1:]) 
