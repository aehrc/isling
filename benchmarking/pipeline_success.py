#### COMPARES THE RESULTS OF THE PIPELINE WITH THE KNOWN INTEGRATIONS #### 

###import libraries 
import pandas as pd 
import argparse
import sys
import numpy as np

###main 
def main(argv): 

	#get arguments 
	parser = argparse.ArgumentParser(description='compares viral integrations identified by pipeline with known viral integrations')
	parser.add_argument('--pipeline', help='output file from pipeline', required = True) 
	parser.add_argument('--ints', help='integrations applied file', required = True)
	parser.add_argument('--read_ints', help='file with files containing information on reads containing viral DNA', required = False) #TODO change this to a requirement 
	args = parser.parse_args()  
	
	#read in integration file 
	ints = pd.read_csv(args.ints, header =0, sep='\t') 

	#read in file listing which reads contain viral DNA 
	read_ints = pd.read_csv(args.read_ints, header=0, sep='\t') #TODO uncomment once batch for file has finnished

	#read in file lising the integrations detected by the pipeline 
	pipe_ints = pd.read_csv(args.pipeline, header = 0, sep = '\t')

	true_vreads, pred_vreads = listIDs(read_ints,pipe_ints)
	correct_reads, missed_reads = listSuccess(true_vreads, pred_vreads)
	getStats(pred_vreads,true_vreads,correct_reads)
	
	#filter out ambiguous reads 
	filt_pipe = filterAmbiguous(pipe_ints)

	#try again after filtering 
	true_vreads, pred_vreads = listIDs(read_ints,filt_pipe)
	correct_reads, missed_reads = listSuccess(true_vreads, pred_vreads)
	print("Stats after filtering...")
	getStats(pred_vreads,true_vreads,correct_reads)

	#create a df of the missed reads
	missed_df = readDataFrame(missed_reads, read_ints) #look at whats up with these ones we missed 
	missed_hPos = findMissed(missed_df)
	evaluateMissed(missed_hPos,ints) 

	#print("columns in pipeline file" +pipe_ints.columns.to_series()) 
	#print("columns in reads file" +read_ints.columns.to_series()) 
	#compare the number of entries in each of the files 
	#print("Number of entries in ints file: "+str(len(ints)))
	#print("Number of entries in host_ints file: "+str(len(host_ints)))
	#print("Number of entries in read_ints file: "+str(len(read_ints)))
	#print("Number of entries in pipe_ints file: "+str(len(pipe_ints)))

	
def listIDs(read_ints,pipe_ints): 
	"""Function which makes lists of the IDs of the viral reads and the predicted viral reads""" 
	#create a list of the IDs of the known viral reads 	
	true_vreads = []
	for i in range(len(read_ints)):
		true_vreads.append(read_ints.loc[i]['fragment_id'])  
		
	#create a list of the IDS of the predicted viral reads 
	pred_vreads = []
	for i in range(len(pipe_ints)): 
		pred_vreads.append(pipe_ints.loc[i]['ReadID']) 

	return true_vreads, pred_vreads

def listSuccess(true_vreads, pred_vreads): 
	"""function which creates a list of the IDs successfully predicted by the pipeline and those missed""" 

	correct_reads = [] #list of reads which were accurately predicted 
	missed_reads = [] #list of reads containing viral DNA missed by the pipeline 

	for int in true_vreads: 
		if int in pred_vreads:
			correct_reads.append(int)
		else: 
			missed_reads.append(int) 

	return correct_reads, missed_reads

def getStats(pred_vreads,true_vreads,correct_reads): 
	"""Function which identifies the accuracy which reads containing viral DNA are identified""" 


	totalPred = len(pred_vreads) #number of reads predicted to contain viral DNA 
	allReads = len(true_vreads) #number of reads actually containing viral DNA 

	#%accuracy of reads identified as containing viral DNA actually containing viral DNA
	#"what chance is there that the predicted read is actually viral"
	predAcc = (len(correct_reads)/totalPred)*100
	print("% OF PIPELINE VIRAL READS CONTAINING VIRAL DNA: "+str(int(predAcc))+"%")

	#% of integrations which were correctly predicted
	#"How good were we at making predictions"  
	predHit = (len(correct_reads)/allReads)*100
	print("% OF VIRAL READS IDENTIFIED BY PIPELINE: "+str(int(predHit))+"%")

	#% of predicitions which were incorrect 
	#"What chance is there that the predicted read is not viral"  
	falsePred = ((totalPred-len(correct_reads))/totalPred)*100
	print("% OF READS FALSELY IDENTIFIED TO CONTAIN VIRAL DNA: "+str(int(falsePred))+"%")  

	#% predictions missed by pipline
	#"What percentage of integrated reads does our pipeline miss"   
	#predMiss = (len(missed_reads)/allReads) 
	predMiss = ((allReads-len(correct_reads))/allReads)*100 
	print("% OF READS READS MISSED BY THE PIPELINE: "+str(int(predMiss))+"%")  
	return predAcc, predHit, falsePred, predMiss 


def filterAmbiguous(pipe_ints): 
	"""Filters out integrated reads which are ambiguous for the host or viral sequence""" 
	#list of the indexes of ambiguous reads 
	ambiguous_idx = []

	for i in range(len(pipe_ints)): 
		if pipe_ints['HostPossibleAmbiguous'][i] == 'yes' or pipe_ints['ViralPossibleAmbiguous'][i] == 'yes': 
			ambiguous_idx.append(i) 
			
	#drop the ambiguous rows 	
	filt_pipe = pipe_ints.drop(pipe_ints.index[ambiguous_idx])

	#reindex fil_pipe for ease of use later 
	filt_pipe = filt_pipe.reset_index(drop=True) 

	#report the % remaining after filtering 
	filt_per = (len(filt_pipe)/len(pipe_ints))*100 
	print(str(int(filt_per))+"% of integrations remaining after filtering reads ambiguous for host or vector") 
 
	return filt_pipe

def readDataFrame(read_list, read_ints): 
	"""Creates a dataframe of a subset of the missed reads using a list of read IDs""" 
	#create a list of the indexes to be included in the new dataframe 
	new_idx = []
	
	read_ints = read_ints.set_index('fragment_id') 
	
	#create df with only the entries in read_list
	new_df = read_ints.loc[read_list]

	#reset index on the new dataframe 
	new_df['fragment_id'] = new_df.index 
	new_df = new_df.reset_index(drop=True)
	del new_df['Unnamed: 0']

	return new_df

def findMissed(missed_df): 
	"""Looks up the integrations missed by the pipeline enabling subsequent determination of what kind of reads we missed""" 

	missed_column = missed_df['hPos'] 

	#list of integration missed by pipeline
	missed_hPos = []
 
	for i in range(len(missed_column)): 	
		values = str(missed_column[i][1:-1])
		values = values.split(", ")
		values = list(map(int,values))
		for value in values:
			#ensure that the same hPos value is not added more than once  
			if value not in missed_hPos:
				missed_hPos.append(value) 
	return missed_hPos 
	 
def evaluateMissed(missed_hPos,ints):
	"""Attempt to find out why the integrations were missed by the pipeline"""
	#drop the columns from ints we do not care about to make it less cumbersome
	ints = ints.drop(['hChr', 'virus'], axis=1) 

	#create a dictionary of properities of each of the integrations
	hPos_dict = ints.set_index('hPos').to_dict()
 
	#check out the different properities to see if we can find any particular issues 
	rearrange = [hPos_dict.get('rearrangement').get(missed) for missed in missed_hPos]
	deletion = [hPos_dict.get('deletion').get(missed) for missed in missed_hPos]
	num_fragments = [hPos_dict.get('num_fragments').get(missed) for missed in missed_hPos]

	#get statistics on the different properties 
	assessRearrange(rearrange) 
	assessDeletion(deletion)
	assessFragments(num_fragments)  

	#TODO modify read_integrations code to say what type of overlap is present in each read 
	left_junction = [hPos_dict.get('leftJunction').get(missed) for missed in missed_hPos]
	right_junction = [hPos_dict.get('rightJunction').get(missed) for missed in missed_hPos]

	#TODO look at the length of the integration (short integrations may be more likely to be missed)
	#actual length of integration or just the part in the read  
 	
def assessRearrange(rearrange): 
	"""Tells us how likely it is that rearrangees are the cause of misses""" 

	#find the number of missed integrations which have been rearranged 
	rearr_true = rearrange.count(True) 

	missed_rearr = (rearr_true/len(rearrange))*100 
	
	print(str(missed_rearr)+"% of the integrations missed were rearraged")   

def assessDeletion(deletion): 
	"""Tells us how likely it is that deletions ar ethe cause of misses""" 
	
	#find the number of missed integrations which had deletions
	del_true = deletion.count(True) 

	missed_rearr = (del_true/len(deletion))*100

	print(str(del)+"% of the integrations missed were rearranged")

def assessFragments(num_fragments): 
	"""Tells us what number of fragments the missed integrations were divided into""" 
	min_frag = min(num_fragments) 
	max_frag = max(num_fragments)


	for i in range(min_frag, max_frag+1):
		count = num_fragments.count(i)
		frag_freq = (count/len(num_fragments))*100 
		print(str(frag_freq)+"% of missed integrations were broken into "+str(i)+" fragments"

if __name__ == "__main__":
	main(sys.argv[1:])
