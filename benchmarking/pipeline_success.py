#### COMPARES THE RESULTS OF THE PIPELINE WITH THE KNOWN INTEGRATIONS #### 

###import libraries 
import pandas as pd 
import argparse
import sys
import numpy as np
from pylab import savefig
import matplotlib as plt
import seaborn as sns 
sns.set
plt.use('cairo') 

###main 
def main(argv): 

	#get arguments 
	parser = argparse.ArgumentParser(description='compares viral integrations identified by pipeline with known viral integrations')
	parser.add_argument('--pipeline', help='output file from pipeline', required = True) 
	parser.add_argument('--ints', help='integrations applied file', required = True)
	parser.add_argument('--read_ints', help='file with files containing information on reads containing viral DNA', required = True)  
	args = parser.parse_args()  
	
	print("Starting") 
	#read in integration file 
	ints = pd.read_csv(args.ints, header =0, sep='\t')
	#filter out the false integrations  
	ints = filterFalse(ints)

	#read in file listing which reads contain viral DNA 
	read_ints = pd.read_csv(args.read_ints, header=0, sep='\t') 
	#filter our file to leave only chimeric reads 
	read_ints = filterChimeric(read_ints)

	#filter our file to remove short reads 
	read_ints = filterLength(read_ints,20)  

	#read in file lising the integrations detected by the pipeline 
	pipe_ints = pd.read_csv(args.pipeline, header = 0, sep = '\t')

	#filter out ambiguous reads 
	filt_pipe = filterAmbiguous(pipe_ints)
 
	true_vreads, pred_vreads = listIDs(read_ints,filt_pipe)
	correct_reads, missed_reads = listSuccess(true_vreads, pred_vreads)
	print("Stats after filtering...")
	getStats(pred_vreads,true_vreads,correct_reads)

	#create a df of the missed reads
	missed_df = readDataFrame(missed_reads,read_ints) #look at whats up with these ones we missed 
	missed_hPos = findMissed(read_ints, correct_reads)
	evaluateMissed(missed_hPos,ints)

	#analyse length of missed reads 
	assessLength(read_ints, missed_reads)

	#get the reads that we correctly identified so that we can look for commonalities
	#missed_df = readDataFrame(correct_reads, read_ints) 
	 
	assessMissedInts(missed_hPos, ints)

	#get the actual types of junctions 
	#assessOverlap(read_ints, missed_reads)  

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
		if read_ints.loc[i]['fragment_id'] not in true_vreads:  #TODO check if makes results better else remove 
			true_vreads.append(read_ints.loc[i]['fragment_id'])  
		
	#create a list of the IDS of the predicted viral reads 
	pred_vreads = []
	for i in range(len(pipe_ints)):
		if pipe_ints.loc[i]['ReadID'] not in pred_vreads: #TODO check otherwise remove 
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

def filterChimeric(read_ints): 
	"""Process of filtering the entire dataframe was inefficient, aim to filter differently""" 
	#list of the indexes of the reads which are not chimeric
	nonchimeric_idx = []

	for i in range(len(read_ints)): 
		if read_ints['left_read'][i] == 'viral' and read_ints['right_read'][i] == 'viral': 
			nonchimeric_idx.append(i) 
			
	#drop the nonchimeric rows 
	filt_reads = read_ints.drop(read_ints.index[nonchimeric_idx])

	#reindex the filtered reads
	filt_reads = filt_reads.reset_index(drop=True) 

	return filt_reads


def filterLength(read_ints, min_len): 
	"""Function which filters the dataframe to contain only reads with more than a set amount of viral DNA. Pipeline can only detect viral DNA if there is more than 20 bp of viral DNA""" 

	#list of the reads which have less than the min_len of base pairs 
	short_idx = []
	
	for i in range(len(read_ints)):
		#take the max of the reads of the fragment  
		v_len = max(read_ints['left_read_amount'][i], read_ints['right_read_amount'][i])
		if v_len <= min_len: 
			short_idx.append(i) 

	#drop the short rows 
	filt_reads = read_ints.drop(read_ints.index[short_idx]) 
	
	#redinex the filtered reads 
	filt_reads = filt_reads.reset_index(drop=True)

	#report the filtering  
	rem = (len(filt_reads)/len(read_ints))*100
	print("After filtering out reads with less than " +str(min_len) +",\n base pair of viral DNA, "+ str(rem)+"% of reads remain.") 

	return filt_reads

	#report the amount of reads remaining after filtering 

def filterFalse(ints): 
	"""Remove integrations which were unsuccessful""" 
	
	#list of the unsuccessful integrations 
	false_idx = []
	
	for i in range(len(ints)): 
		if ints['inserted'][i] == False: 
			false_idx.append(i)

	#drop the false rows 
	filt_ints = ints.drop(ints.index[false_idx]) 

	#reindex the filtered reads
	filt_ints = filt_ints.reset_index(drop=True) 

	return filt_ints 

def readDataFrame(read_list, read_ints): 
	"""Creates a dataframe of a subset of  reads using a list of read IDs""" 
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


def findMissed(read_ints, correct_reads): 


	"""Looks up the integrations missed by the pipeline enabling subsequent determination of what kind of reads we missed""" 

	#find the maximum number of integrations that could be found from the reads if the pipeline was perfect 
	detectable_hPos = []
		
	for i in range(len(read_ints)):
		values = str(read_ints['hPos'][i][1:-1])
		values = values.split(", ")
		values = list(map(int,values))
		for value in values:
			if value not in detectable_hPos:  
				detectable_hPos.append(value) 
	print('Max number of integrations detectable from reads: '+str(len(detectable_hPos))) 
	#TODO look at which ones were not detected - potentially can be improved by increasing the coverage

	#look at how many of these integrations we cover with our results from the pipeline 
	#make a dataframe of the correct reads 
	correct_df = readDataFrame(correct_reads, read_ints) 

	#find the number of integrations that the pipeline detected 
	pipeline_hPos = []
	
	for i in range(len(correct_df)): 
		values = str(correct_df['hPos'][i][1:-1])
		values = values.split(", ")
		values = list(map(int,values))
		for value in values:
			if value not in pipeline_hPos: 
				pipeline_hPos.append(value) 

	print("Number of integrations detectable from the pipeline: "+str(len(pipeline_hPos))) 

	acc = (len(pipeline_hPos)/len(detectable_hPos))*100

	print("Of the integrations which were present in the reads, "+str(acc)+"% could be detected") 
   
	#create a column of the missed integrations so we can evaluate them
	missed_hPos =  []
	for hPos in detectable_hPos: 
		if hPos not in pipeline_hPos: 
			missed_hPos.append(hPos) 


	return missed_hPos 
	 
def evaluateMissed(missed_hPos,ints):
	"""Attempt to find out why the integrations were missed by the pipeline""" 
	#create a dictionary of properities of each of the integrations
	hPos_dict = ints.set_index('hPos').to_dict()
 
	#check out the different properities to see if we can find any particular issues 
	rearrange = [hPos_dict.get('rearrangement').get(missed) for missed in missed_hPos]
	deletion = [hPos_dict.get('deletion').get(missed) for missed in missed_hPos]
	num_fragments = [hPos_dict.get('num_fragments').get(missed) for missed in missed_hPos]

	#get statistics on the different properties 
	assessRearrange(rearrange) 
	assessDeletion(deletion)
	assessFragments(num_fragments)  #TODO compare these to total amount which were deletions etc. 

	#report what percentage of the integrations we found 
	missed_int = (len(missed_hPos)/len(ints))*100 
	print('\nThe percentage of integrations missed was '+str(missed_int)+"%") 

def assessMissedInts(missed_hPos, ints): 
	"""Look at the properities of the integrations which were not detected by the pipeline""" 
	#make dictionary of integrations  
	hPos_dict = ints.set_index('hPos').to_dict()
	hPos_keys = list(hPos_dict.get('vBases').keys())

	#get a list of the lengths of the integrations 
	all_lengths = []

	for i in range(len(hPos_keys)): 
		seq = hPos_dict.get('vBases')[hPos_keys[i]]  
		length = len(seq) 
		all_lengths.append(length)
 
	hPos_length = dict(zip(ints['hPos'],all_lengths)) 
	
	#save distribution of all integration lengths  #TODO consider the number of fragments here 
	ax = sns.distplot(all_lengths, bins=20)
	ax.set_title("Histogram of the length of all integrations") 
	ax.set(xlabel='length of integration', ylabel = 'Density')  
	fig = ax.get_figure()
	fig.savefig("integration_lengths.pdf")
	print("Distribution of all integrations saved to integration_lengths.pdf") 

	#save distribution of the missed integration lengths
	missed_ints = [hPos_length.get(missed_int) for missed_int in missed_hPos] 
	ax = sns.distplot(missed_ints, bins=20)
	ax.set_title("Histogram of the length of missed integrations") 
	ax.set(xlabel='length of integration', ylabel = 'Density')  
	fig = ax.get_figure()
	fig.savefig("missed_integration_lengths.pdf")
	print("Distribution of all integrations saved to missed_integration_lengths.pdf")
	
def assessRearrange(rearrange): 
	"""Tells us how likely it is that rearrangees are the cause of misses""" 

	#find the number of missed integrations which have been rearranged 
	rearr_true = rearrange.count(True) 

	missed_rearr = (rearr_true/len(rearrange))*100 
	
	print("\n" + str(missed_rearr)+"% of the integrations missed were rearraged\n")   

def assessDeletion(deletion): 
	"""Tells us how likely it is that deletions ar ethe cause of misses""" 
	
	#find the number of missed integrations which had deletions
	del_true = deletion.count(True) 

	missed_del = (del_true/len(deletion))*100

	print(str(missed_del)+"% of the integrations missed had deletions\n")

def assessFragments(num_fragments): 
	"""Tells us what number of fragments the missed integrations were divided into""" 
	min_frag = min(num_fragments) 
	max_frag = max(num_fragments)


	for i in range(min_frag, max_frag+1):
		count = num_fragments.count(i)
		frag_freq = (count/len(num_fragments))*100 
		print(str(frag_freq)+"% of missed integrations were broken into "+str(i)+" fragments")


def getOverlap(read_ints):
	"""Tells us what type of junctions the filter/filtered reads had (from the true file not the pipeline)
Inormation corresponds to read ID""" 
	#get the requried columns 
	left_read= read_ints['left_read'].values
	right_read = read_ints['right_read'].values
	left_junc = read_ints['left_junc'].values
	right_junc = read_ints['right_junc'].values
	int_ID = read_ints['fragment_id'].values

	#column which says the type of overlaps 
	overlap_type = []

	for i in range(len(read_ints)):
		if left_read[i] == "viral" and right_read[i] == 'host' or right_read[i] == 'viral' and left_read[i] == 'host': 
			overlap_type.append("discordant")
		elif left_read[i] == 'chimeric' and right_read[i] == 'host': 
			overlap_type.append(right_junc[i])
		elif left_read[i] == 'chimeric' and right_read[i] == 'viral': 
			overlap_type.append(left_junc[i]) 
		elif left_read[i] == 'host' and right_read[i] == 'chimeric':
			overlap_type.append(left_junc[i])
		elif left_read[i] == 'viral' and right_read[i] == 'chimeric': 
			overlap_type.append(right_junc[i])
		elif left_read[i] == 'short' or right_read[i] == 'short': 
			overlap_type.append("SHORT READ") #TODO think about whatto do with this 

	#make into a dictionary - easy to make comparison 
	overlap_dict = dict(zip(int_ID, overlap_type)) 

	return overlap_dict 
	
def assessLength(read_ints, missed_reads): 
	"""Compare the number of viral base pairs in the missed reads and in alfg
l reads"""
	#get the required columns
	int_ID = read_ints['fragment_id'].values
	first_len = read_ints['left_read_amount'].values
	second_len = read_ints['right_read_amount'].values

	#find the length of virus in each read - excluding reads consisting of whole virus and no virus 
	viral_len = [max(first_len[i],second_len[i]) for i in range(len(read_ints))]	
	
	#create a dictionay
	len_dict = dict(zip(int_ID, viral_len)) 
	#get the lengths for missed reads 
	missed_len = [len_dict.get(missed) for missed in missed_reads] 
	
	#save distibution of length of the missed data 
	ax = sns.distplot(missed_len, hist=False)
	ax.set_title("Distrubtion of the length of missed integrations") 
	ax.set(xlabel='length of missed integration', ylabel = 'Density')  
	fig = ax.get_figure()
	fig.savefig("missedread_len_distplot.pdf")
	print("Distribution of the length of the missed reads saved to missedread_len_distplot.pdf") 

	
#def assessIncorrectPredictions #TODO look at the integrations which were falsely detected by the pipeline 

def assessOverlap(read_ints, missed_reads): 
	"""Looks at what type of overlap the reads missed by the pipeline have. Takes a dataframe of all of the reads and a list of the IDs of the missed reads""" 
	#dictionary of overlap types for each read 
	overlap_dict = getOverlap(read_ints)
	
	#column of the overlap of the missed_reads
	overlap = [overlap_dict.get(missed) for missed in missed_reads] 
	
	#find the instances of the gap types 
	gap = overlap.count('gap')
	overlap_ = overlap.count('overlap')
	none = overlap.count('none') #clean 
	other = overlap.count('nan') 

	#print the instances of each 
	gapFreq = (gap/len(missed_reads))*100 
	print(str(gapFreq)+"\n% of missed reads had a gap") 
	overlapFreq = (overlap_/len(missed_reads))*100
	print(str(overlapFreq)+"% of missed reads had an overlap") 
	cleanFreq = (none/len(missed_reads))*100
	print(str(cleanFreq)+"% of missed reads had clean junctions") 
	otherFreq =(other/len(missed_reads))*100
	print(str(otherFreq)+"%of missed reads had a junction other than those above") 
	

	

def compareOverlap(read_ints, pipe_ints): 
	"""Compares whether the pipeline has correctly predicted the gap type of the read. Takes list of predicted"""

	#get the IDs of the pipeline reads
	pipe_reads = pipe_ints['ReadID'].values 

	#get the overlap type of the pipeline reads
	pipe_overlap = pip_ints['Type'].values
 
	#create a dictionary of the overlap types for each read 
	overlap_dict = getOverlap(read_ints)
	
	#get the known overlap type 
	actual_overlap = [overlap_dict.get(pipeline) for pipeline in pipe_reads] 
	#TODO handle reads in pipeline but don't actually contain viral DNA 

	#make comparisons between pipeline overlaps and known overlaps 
	correct = 0 #counter for the number of correct estimations 
	for i in range(len(pipe_reads)): 
		if pipe_overlaps[i] == actual_overlap[i]: 
			correct = correct +1 

	#caluclate the accuracy of the pipeline to find the type of junction
	accuracy = (correct/len(pipe_reads))*100 
	print("Pipeline predicted the type of junction correctly "+str(accuracy)+"% of the time:") 

	#we also want to know which ones we got wrong 
	#find the distribtuion of rearrangements and deletions so that these can be compared - plot using seaborn? 


if __name__ == "__main__":
	main(sys.argv[1:])
