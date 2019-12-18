#### COMPARES THE RESULTS OF THE PIPELINE WITH THE KNOWN INTEGRATIONS #### 

###import libraries 
import pandas as pd 
import argparse
import sys

###main 
def main(argv): 

	#get arguments 
	parser = argparse.ArgumentParser(description='compares viral integrations identified by pipeline with known viral integrations')
	parser.add_argument('--pipeline', help='output file from pipeline', required = True) 
	parser.add_argument('--ints', help='integrations applied file', required = True)
	parser.add_argument('--host_ints', help='file listing regions with viral DNA', required = True)
	parser.add_argument('--read_ints', help='file with files containing information on reads containing viral DNA', required = False) #TODO change this to a requirement 
	args = parser.parse_args()  
	
	#read in integration file 
	ints = pd.read_csv(args.ints, header =0, sep='\t') 

	#read in file listing regions with viral DNA 
	host_ints = pd.read_csv(args.host_ints,header=0,sep='\t')

	#read in file listing which reads contain viral DNA 
	read_ints = pd.read_csv(args.read_ints, header=0, sep='\t') #TODO uncomment once batch for file has finnished

	#read in file lising the integrations detected by the pipeline 
	pipe_ints = pd.read_csv(args.pipeline, header = 0, sep = '\t')

	true_vreads, pred_vreads = listIDs(read_ints,pipe_ints)
	correct_reads, missed_reads = listSuccess(true_vreads, pred_vreads)
	getStats(pred_vreads,true_vreads,correct_reads)
	
	
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
	predMiss = ((allReads-len(correct_reads))/allReads) #compare TODO
	print("% OF READS READS MISSED BY THE PIPELINE: "+str(int(predMiss))+"%")  

	return predAcc, predHit, falsePred, predMiss 

def findReadType(): 
	"""Function which identifies the accuracy which the type of read is identified by the pipeline""" 
		

def findStartStop(): 
	"""Function which identifies the accuracy which the start and stop points are identified by the pipeline""" 





if __name__ == "__main__":
	main(sys.argv[1:])
