### finds the seq, CIGAR, XA and SA of reads ##

## pipeline appears to produce more false positives and false negatives for reads containing integrations with overlaps 
## to determine whether issues are a result of error in the intial insert_virus.py or read_integrations.py scripts can can further 
## analyse these false reads using the corresponding seq, CIGAR, XA (secondary alignment) and SA (supplementary alignment) 

##import libraries 
import pandas as pd
from simplesam import Reader
import argparse


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

	#get the bam file


	#before continuing add code to pipeline_success to obtain the false positive and false negative IDs 


if __name__ == "__main__":
	main(sys.argv[1:]) 





