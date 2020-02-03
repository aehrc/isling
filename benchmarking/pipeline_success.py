#### COMPARES THE RESULTS OF THE PIPELINE WITH THE KNOWN INTEGRATIONS #### 

###import libraries 
import pandas as pd 
import argparse
import sys
import numpy as np
from pylab import savefig
import matplotlib as plt
import seaborn as sns
import os   
sns.set
plt.use('cairo') 

###main 
def main(argv): 

	#get arguments 
	parser = argparse.ArgumentParser(description='compares viral integrations identified by pipeline with known viral integrations')
	parser.add_argument('--pipeline', help='output file from pipeline', required = True) 
	parser.add_argument('--ints', help='integrations applied file', required = True)
	parser.add_argument('--viral_reads', help='file with reads containing information on reads containing viral DNA', required = True)
	parser.add_argument('--all_reads', help = 'file containing all reads IDs', required = True) 
	args = parser.parse_args()  
	
	
	print("Output to be saved to /evaluate_pipeline_output/output.txt") 

	#create directory to save output
	os.makedirs('evaluate_pipeline_output', exist_ok = True)

	#save output in terminal to a file 
	f = open( "evaluate_pipeline_output/output.txt", "w")
	sys.stdout = f  
	
	print("STARTING...", flush = True) 
	
	#read in integration file 
	ints = pd.read_csv(args.ints, header =0, sep='\t')
	#filter out the false integrations - we don't care about this for detecting integrations  
	ints = filterFalse(ints)

	#read in file listing all reads 
	all_reads = pd.read_csv(args.all_reads,header = 0, sep = '\t')
	all_IDs = list(set(all_reads["fragment_id"]))
	all_IDs = [i.replace('chr', '') for i in all_IDs] 	


	#read in file listing which reads contain viral DNA 
	viral_reads = pd.read_csv(args.viral_reads, header=0, sep='\t')
	viral_reads['fragment_id'] = [x.replace("chr","") for x in viral_reads['fragment_id']]
	#pipeline can only detect a read as being viral if there is at least 20 bases each of host and virus DNA
	all_reads = filterLength(viral_reads, 20)  

	#read in file listing the integrations detected by the pipeline 
	pipe_ints = pd.read_csv(args.pipeline, header = 0, sep = '\t')
	print("Number of reads detected by pipeline: " +str(len(pipe_ints)), flush = True)

	#look at the different types of filtering 
	compareFilters(pipe_ints,all_reads, all_IDs)

	#look for what ratio of integrations we captured with the detected reads
	"""
	actual_Vreads, pred_Vreads, actual_NVreads, pred_NVreads = listIDs(all_reads, pipe_ints, all_IDs)
	detected_Vreads, undetected_Vreads, detected_NVreads, undetected_NVreads = listSuccess(actual_Vreads, actual_NVreads, pred_Vreads, pred_NVreads, False) 
	missed_hPos = findMissed(all_reads, detected_Vreads)
	"""

	f.close()

def filterList(pipe_ints, ID_list):
	"""Filters the set of pipeline reads so that only mapped reads are included. Created to look at whether pipeline is including unmapped reads resulting in high false positive rate""" 

	#create list of indexes to be filtered
	filter_idx = []

	pipe_IDs = list(set(pipe_ints["ReadID"]))
	ID_list = [i.replace("chr"," ") for i in ID_list]
	ID_list = list(set(ID_list))

	print(pipe_IDs[1:5], flush = True) 
	print(ID_list[1:5], flush = True) 
	
	for i in range (len(pipe_IDs)): 
		if pipe_IDs[i] in ID_list: 
			filter_idx.append(i)

	#drop the unmapped rows 
	filt_pipe = pipe_ints.drop(pipe_ints.index[filter_idx])
	
	#reindex the filtered pipeline
	filt_pipe = filt_pipe.reset_index(drop=True) 

	print("Entries before filtering: "+str(len(pipe_ints)), flush = True)
	print("Entries after filtering: "+str(len(filt_pipe)), flush = True)
	return filt_pipe
			

def listIDs(viral_reads, pipe_ints, all_IDs): 
	"""Create lists of the predicted and actual viral and non-viral reads""" 

	#List all IDs 
	all_reads = list(set(all_IDs))
	print("Number of reads: "+str(len(all_reads)), flush = True)

	#list pipeline predictions
	pipe_IDs =list(set(pipe_ints["ReadID"]))
	#remove space at the start of the ID
	pipe_IDs = [i.replace(" ", "") for i in pipe_IDs]  
	
	#List actual viral reads 
	actual_Vreads = list(set(viral_reads['fragment_id']))
	actual_Vreads = [i.replace("chr","") for i in actual_Vreads] 
	print("Number of viral reads: "+str(len(actual_Vreads)), flush = True) 

	#List viral reads predicted by the pipeline 
	pred_Vreads = list(set(pipe_IDs).intersection(set(all_reads)))
	print("Number of viral reads predicted by the pipeline: "+str(len(pred_Vreads)), flush = True)

	#List viral reads appearing in pipeline which are unknown 
	pred_unk = list(set(pipe_IDs)-set(all_reads)) 
	print("Reads filtered in processing (ie low quality): "+str(len(pred_unk)), flush = True) 

	#List actual non-viral reads 
	actual_NVreads = list(set(all_reads) - set(actual_Vreads)) 
	print("Number of non-viral reads: "+str(len(actual_NVreads)), flush = True) 

	#list non-viral reads predicted by the pipeline 
	pred_NVreads = list(set(all_reads)-set(pred_Vreads))
	print("Number of predicted non-viral reads: " +str(len(pred_NVreads)), flush = True) 

	print("Sum of actual reads: "+str(len(actual_Vreads)+len(actual_NVreads)),flush = True)
	print("Sum of predicted reads: "+str(len(pred_Vreads)+len(pred_NVreads)),flush = True)

	return actual_Vreads, pred_Vreads, actual_NVreads, pred_NVreads 


def listSuccess(actual_Vreads, actual_NVreads, pred_Vreads, pred_NVreads, save): 
	"""function which creates a list of the IDs successfully predicted by the pipeline and those missed. actual Vreads is a list of the reads known to be viral and pred_Vreads is a list of the reads predicted by the pipeline to contain viral DNA. Can save false positives and negatives to file if save == True.""" 

	#find how many of the viral reads were/were not detected by the pipeline
	detected_Vreads = list(set(actual_Vreads).intersection(pred_Vreads))
	print("\nTrue Positives: "+str(len(detected_Vreads)))
	undetected_Vreads = list(set(actual_Vreads).intersection(pred_NVreads))
	print("False Negatives: "+str(len(undetected_Vreads))) 

	#find how many of the non-viral reads were/were not detected by the pipeline 
	detected_NVreads = list(set(actual_NVreads).intersection(pred_Vreads))
	print("False Positives: "+str(len(detected_NVreads))) 
	undetected_NVreads = list(set(actual_NVreads).intersection(pred_NVreads))
	print("True Negatives: "+str(len(undetected_NVreads)))

	pred_pos = len(detected_Vreads)+len(detected_NVreads)
	pred_neg = len(undetected_Vreads)+len(undetected_NVreads)

	print("Total number of positive predictions: " + str(pred_pos)) #debugging
	print("Total number of negative predictions: "+str(pred_neg)) #debugging 
	
	if save == True: 
		#save the false positive reads to file
		with open('evaluate_pipeline_output/false_positive_IDs.txt', 'w') as f: 
			for item in detected_NVreads: 
				f.write("%s\n" % item)
		#save the false negative reads to file 
		with open('evaluate_pipeline_output/false_negative_IDs.txt', 'w') as f: 
			for item in undetected_Vreads: 
				f.write("%s\n" % item)
		print("False positive reads saved to 'evaluate_pipeline_output/false_positive_IDs.txt'", flush = True) 
		print("False negative reads saved to 'evaluate_pipeline_output/false_negative_IDs.txt'", flush = True)

	return detected_Vreads, undetected_Vreads, detected_NVreads, undetected_NVreads

def findStats(detected_Vreads, undetected_Vreads, detected_NVreads, undetected_NVreads,save): 
	"""Get statisitics on the the success of the pipeline""" 

	# find the number of TP, FN, FP and TN 
	TP = len(detected_Vreads)
	FN = len(undetected_Vreads)
	FP = len(detected_NVreads)
	TN = len(undetected_NVreads)
	
  	#find corresponding rates
	TPR = TP/(TP+FP)*100
	TNR = TN/(TN+FN)*100
	FPR = FP/(TP+FP)*100 
	FNR = FN/(TN+FN)*100 

	#construct confusion matrix
	conf_mat = np.array([[TPR,FPR],[FNR, TNR]])
	consfusionMatrix(conf_mat)

	#accuracy - chance of making a correct prediction
	acc = ((TP+TN)/(TP+TN+FP+FN))*100
	print("Prediction accuracy: "+str(acc)+"%",flush=True)

	#precision - chance that positive test result is actually a positive test result 
	prec = (TP/(TP+FP))*100
	print("Prediction precision: "+str(prec)+"%",flush = True) 

	#sensitivity 
	sens = (TP/(TP+FN))*100
	print("Prediction sensitivity: "+str(sens) +"%", flush = True)

	#specificity 
	spec = (TN/(TN+FP))*100
	print("Prediction Specificity: "+str(spec)+"%", flush = True)

	#F-measure
	F1 = ((2*prec*sens)/(prec+sens))
	print("F1 measure: "+str(F1)+"%",flush = True) 

	#create a dataframe of the values and then these can be appended to create a graph
	stats = pd.DataFrame({"Accuracy": [acc], "Precision": [prec],"Sensitivity": [sens], "F measure": [F1], "Specificity": [spec]}) 

	#save the confusion matrix. Useful if want to make a separate notebook later to look at results 
	conf_df = pd.DataFrame({"TP":[TP], "TN": [TN], "FP": [FP], "FN": [FN]})

	return stats, conf_df

def consfusionMatrix(conf_mat): 
	"""Makes and saves a confusion matrix of the accuracy of the pipeline""" 
	
	fig, ax = plt.pyplot.subplots() 
	axis_labels = ["Viral", "Non-viral"]
	ax = sns.heatmap(conf_mat, annot= conf_mat, fmt = '.1f',xticklabels = axis_labels, 
yticklabels = axis_labels, annot_kws = {"size":15}, cmap = 'Blues')
	for t in ax.texts: t.set_text(t.get_text()+"%")
	ax.set_xlabel('Actual',fontsize = 15)
	ax.set_ylabel('Pipeline Predicted',fontsize = 15)
	fig = ax.get_figure()
	fig.savefig("evaluate_pipeline_output/confusion_matrix.pdf")
	print("Confusion matrix saved at evaluate_pipeline_output/confusion_matrix.pdf", flush = True)

def barGraph(stats, xlabel):
	"""Creates a bar graph of the accuracy, sensitivity and precision of the pipeline. Append all of the stats into one dataframe before parsing."""

	plt.style.use("ggplot")
	
	#find the number of groups to graph
	n = len(stats) 
	
	fig, ax = plt.pyplot.subplots()
	plt.pyplot.gcf().subplots_adjust(bottom = 0.5) 
	index = np.arange(n) 
	bar_width = 0.2
	opacity = 0.9
	ax.bar(index, stats["Accuracy"], bar_width, alpha = opacity, color = 'r', label = "Accuracy")
	ax.bar(index+bar_width, stats["Precision"], bar_width, alpha = opacity, color = 'b', label = "Precision")
	ax.bar(index+2*bar_width, stats["Sensitivity"], bar_width, alpha = opacity, color = 'g', label = "Sensitivity")
	ax.bar(index+3*bar_width, stats["F measure"], bar_width, alpha = opacity, color = 'darkorange', label = "F measure")
	ax.set_xlabel(xlabel)
	ax.set_xticks(index+bar_width) 
	ax.set_xticklabels(stats.index.values.astype(str), rotation = 45, ha = 'right') #make sure the index of the stats dataframe is named
	ax.legend(ncol = 4)
	fig = ax.get_figure()
	fig.savefig("evaluate_pipeline_output/bar_graph_tests.pdf") 
	print("Bar graph of tests saved at evaluate_pipeline_output/confusion_matrix.pdf", flush = True) 

	
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
	print("{:.2f}% of integrations remaining after filtering reads ambiguous for host or vector".format(filt_per), flush = True) 
 
	return filt_pipe

def currentFiltering(pipe_ints): 
	"""Simulates the current filtering which Suzanne is using to identify which reads are real""" 
	
	#list of the indexes to filter 
	filter_idx = []

	#max edit distance - can be decreased to be more conservative 
	edit_dist =  5 

	#max ambiguous bases - can be decreased to be more conservative 
	max_ambig = 20 

	#intialise counters for each filtering 
	hed = 0 #host edit distance 
	ved = 0 #viral edit distance
	nab = 0 #noambiguousambases 
	otd = 0 #overlap type = discordant
	pvr = 0 #possible vector rearrangement
	pht = 0 #possible host translocation
	
	
	for i in range(len(pipe_ints)): 
		# identify reads with high HostEditDistance 
		if pipe_ints['HostEditDist'][i] > edit_dist: 
			filter_idx.append(i)
			hed += 1
		# identify reads with high ViralEditDistance 
		if pipe_ints['ViralEditDist'][i] > edit_dist: 
			filter_idx.append(i)
			ved += 1
		#filter if too many ambiguous bases
		if pipe_ints['NoAmbiguousBases'][i] == "?" or int(pipe_ints['NoAmbiguousBases'][i]) > max_ambig:
			filter_idx.append(i)
			nab += 1 
		#filter discordant reads 
		if pipe_ints['OverlapType'][i] == 'discordant': #TODO unsure of this 
			filter_idx.append(i) 
			otd += 1
		#filter vector rearrangement 
		if pipe_ints['PossibleVectorRearrangement'][i] == 'yes':
			filter_idx.append(i)
			pvr += 1
		#filter host rearrangement 
		if pipe_ints['PossibleHostTranslocation'][i] == 'yes': 
			filter_idx.append(i)
			pht += 1
		

	print(str(hed)+" Reads filtered with HostEditDist <= "+str(edit_dist),flush=True)
	print(str(ved)+" Reads filtered with ViralEditDist <= "+str(edit_dist),flush=True)
	print(str(nab)+" Reads filtered to contain reads only with NoAmbiguousBases <= "+str(max_ambig), flush = True)
	print(str(otd)+" Discordant reads filtered", flush=True)
	print(str(pht)+" Reads filtered with possible host translocations")
	print(str(pvr)+" Reads filtered with possible vector rearrangements")

	#drop the filtered rows 
	filt_pipe = pipe_ints.drop(pipe_ints.index[filter_idx]) 
	
	#reindex filt_pipe for ease of use later 
	filt_pipe = filt_pipe.reset_index(drop=True) 

	#report the % which were filtered out - important for knowing how much data remains 
	filt_per = (len(filt_pipe)/len(pipe_ints))*100 
	print("{:.2f}% of integrations remaining after applying the current filtering".format(filt_per), flush = True) 
 
	return filt_pipe

def compareFilters(pipe_ints,all_reads, all_IDs):
	"""Create function to compare the different types of filters on data and complete statistics after each""" 
	#TODO tidy this up
	
	#Look at host edit distance 
	filter_idx = []
	edit_dist = 5 
	for i in range(len(pipe_ints)): 
		# identify reads with high HostEditDistance 
		if pipe_ints['HostEditDist'][i] > edit_dist: 
			filter_idx.append(i)
	stats, conf_df = filteredStats(filter_idx, pipe_ints, "HostEditDist > 5", all_IDs, all_reads, False)
	#create dataframe to append stats to 
	all_stats = pd.DataFrame(stats)
	all_conf = pd.DataFrame(conf_df) 

	#look at viral edit distance
	filter_idx = []
	edit_dist = 5
	for i in range(len(pipe_ints)): 
		# identify reads with high HostEditDistance 
		if pipe_ints['ViralEditDist'][i] > edit_dist: 
			filter_idx.append(i)
	stats, conf_df = filteredStats(filter_idx, pipe_ints, "ViralEditDist > 5", all_IDs, all_reads, False)
	#add new stats to dataFrame 
	all_stats = all_stats.append(stats, ignore_index = True)
	all_conf = all_conf.append(conf_df, ignore_index = True) 

	#look at total edit distance
	filter_idx = []
	total_dist = 7 
	for i in range(len(pipe_ints)): 
		if pipe_ints['TotalEditDist'][i] > total_dist: 
			filter_idx.append(i)
	stats, conf_df = filteredStats(filter_idx, pipe_ints, "TotalEditDist > 7", all_IDs, all_reads, False)
	#add new stats to dataFrame 
	all_stats = all_stats.append(stats, ignore_index = True)
	all_conf = all_conf.append(conf_df, ignore_index = True) 
	
	#look at noambiguous bases
	filter_idx = []
	max_ambig = 20
	for i in range(len(pipe_ints)):
		if pipe_ints['NoAmbiguousBases'][i] == "?" or int(pipe_ints['NoAmbiguousBases'][i]) > max_ambig:
			filter_idx.append(i)
	stats, conf_df = filteredStats(filter_idx, pipe_ints, "NoAmiguousBases == ? | 20", all_IDs, all_reads, False)
	#add new stats to dataFrame 
	all_stats = all_stats.append(stats, ignore_index = True)
	all_conf = all_conf.append(conf_df, ignore_index = True) 
	
	#look at discordant reads 
	filter_idx = []
	for i in range(len(pipe_ints)):
		if pipe_ints['OverlapType'][i] == 'discordant': 
			filter_idx.append(i)
	stats, conf_df = filteredStats(filter_idx, pipe_ints, "Discordant", all_IDs, all_reads, False)
	#add new stats to dataFrame 
	all_stats = all_stats.append(stats, ignore_index = True)
	all_conf = all_conf.append(conf_df, ignore_index = True)  

	#look at vector rearrangements
	filter_idx = []
	for i in range(len(pipe_ints)):
		if pipe_ints['PossibleVectorRearrangement'][i] == 'yes':
			filter_idx.append(i)
	stats, conf_df = filteredStats(filter_idx, pipe_ints, "PossibleVectorRearrangement", all_IDs, all_reads, False)
	#add new stats to dataFrame 
	all_stats = all_stats.append(stats, ignore_index = False)
	all_conf = all_conf.append(conf_df, ignore_index = False) 

	#look at host rearrangements 
	filter_idx = []
	for i in range(len(pipe_ints)):
		if pipe_ints['PossibleHostTranslocation'][i] == 'yes':
			filter_idx.append(i)
	stats, conf_df = filteredStats(filter_idx, pipe_ints, "PossibleHostTranslocation", all_IDs, all_reads, False)
	#add new stats to dataFrame 
	all_stats = all_stats.append(stats, ignore_index = True)
	all_conf = all_conf.append(conf_df, ignore_index = True) 

	#look at host possible ambiguous 
	filter_idx = []
	for i in range(len(pipe_ints)): 
		if pipe_ints['HostPossibleAmbiguous'][i] == 'yes': 
			filter_idx.append(i) 
	stats, conf_df = filteredStats(filter_idx, pipe_ints, "HostPossibleAmbiguous", all_IDs, all_reads, False) 
	#add new stats to dataFrame 
	all_stats = all_stats.append(stats, ignore_index = True)
	all_conf = all_conf.append(conf_df, ignore_index = True) 

	#look at vector possible ambiguous 
	filter_idx = []
	for i in range(len(pipe_ints)): 
		if pipe_ints['ViralPossibleAmbiguous'][i] == 'yes':
			filter_idx.append(i) 
	stats, conf_df = filteredStats(filter_idx, pipe_ints, 'ViralPossibleAmbiguous', all_IDs, all_reads, False) 
	#add new stats to dataFrame 
	all_stats = all_stats.append(stats, ignore_index = True)
	all_conf = all_conf.append(conf_df, ignore_index = True)

	#look at no filter
	filter_idx = []
	stats, conf_df = filteredStats(filter_idx, pipe_ints, "Nofilter", all_IDs, all_reads, True)
	#add new stats to dataFrame 
	all_stats = all_stats.append(stats, ignore_index = True)
	all_conf = all_conf.append(conf_df, ignore_index = True) 

	#give the dataframes indexes
	idx = ["HostEditDist <= 5", "ViralEditDist <=5", "TotalEditDist <=6" , "NoAmbiguousBases <= 20", "No Discordant", "PossibleVectorRearrangement", "PossibleHostTranslocation", "HostPossibleAmbiguous",'ViralPossibleAmbiguous', "No filter"]
	all_stats.index = idx
	all_conf.index = idx 

	#construct graph
	barGraph(all_stats, "Filtering")

	#save the statistics if we want to analyse them/compare with other results (ie in a notebook) 
	all_stats.to_csv("evaluate_pipeline_output/filtering_stats.csv", sep = '\t', index = True)
	all_conf.to_csv("evaluate_pipeline_output/conf_stats.csv", sep = '\t', index = True)  

def filteredStats(filter_idx, pipe_ints, tag, all_IDs, all_reads, save): 
	"""Makes a dataframe of filtered statistics for a list of indexes to be removed"""
 	
	#remove filtered indexes from pipeline results
	filt_pipe = pipe_ints.drop(pipe_ints.index[filter_idx]) 
	filt_pipe = filt_pipe.reset_index(drop=True)
	print("\n"+str(len(filter_idx))+" "+tag+" filtered" ,flush = True)

	#report statistics of filtering
	actual_Vreads, pred_Vreads, actual_NVreads, pred_NVreads  = listIDs(all_reads,filt_pipe, all_IDs)
	print("Stats after filtering...")
	detected_Vreads, undetected_Vreads, detected_NVreads, undetected_NVreads = listSuccess(actual_Vreads, actual_NVreads, pred_Vreads, pred_NVreads, save)
	stats, conf_df = findStats(detected_Vreads, undetected_Vreads, detected_NVreads, undetected_NVreads, tag)

	return stats, conf_df 

def filterVectorRearrangement(pipe_ints): 
	"""Filter out reads detected by the pipeline which are possible vector rearrangements. This may reduce the false positive rate of the pipeline""" 
	
	#list the indexes of possible vector rearrangements 
	rearr_idx = []

	#obtain indexes of possible vector rearrangements 
	for i in range(len(pipe_ints)):
		if pipe_ints['PossibleVectorRearrangement'][i] == 'yes': 
			rearr_idx.append(i) 

	#drop reads which are possibly vector rearrangements
	filt_pipe = pipe_ints.drop(pipe_ints.index[rearr_idx])
	
	#reindex the df for later use 
	filt_pipe = filt_pipe.reset_index(drop=True) 

	#report the % remaining after filtering 
	filt_per = (len(filt_pipe)/len(pipe_ints))*100 
	print("{:.2f}% of integrations remaining after filtering possible vector rearrangements".format(filt_per), flush= True)
	

	return filt_pipe 

def filterLength(all_reads,min_len):
	"""Function which filters viral reads to contain only reads with a set amount of viral DNA (min_len)""" 

	#set the read length
	read_len = 150 
	
	#make a list of the indexes of the reads we wish to filter
	filt_idx = []

	#loop through and adjust chimeric reads 
	for i in range(len(all_reads)):

		#adjust left read to meet the threshold of 20bp to be chimeric 	
		if all_reads['left_read_amount'][i] > read_len - min_len: 
			all_reads.loc[i,'left_read'] = 'v'
		elif all_reads['left_read_amount'][i] < min_len: 
			all_reads.loc[i,'left_read'] = 'h'

		#adjust right read to meet the threshold of 20bp to be chimeric 
		if all_reads['right_read_amount'][i] > read_len - min_len: 
			all_reads.loc[i,'right_read'] = 'v'
		elif all_reads['right_read_amount'][i] < min_len: 
			all_reads.loc[i,'right_read'] = 'h'

		#if both left and right are viral or host do integration occured here
		if all_reads['right_read'][i] == 'h' and all_reads['left_read'][i] == 'h' or all_reads['right_read'][i] == 'v' and all_reads['left_read'][i] == 'v':
				filt_idx.append(i) 

	#drop the false rows
	filt_reads = all_reads.drop(all_reads.index[filt_idx])
	
	#reindex the filtered reads
	filt_reads = filt_reads.reset_index(drop=True)
	
	#report the filtering 
	#report the filtering  
	rem = (len(filt_reads)/len(all_reads))*100
	print("\nAfter filtering out reads with less than "+str(min_len)+" base pairs of viral DNA, {:.2f} % of reads remain.".format(rem), flush = True)

	return filt_reads
	

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

def readDataFrame(read_list, all_reads): 
	"""Creates a dataframe of a subset of  reads using a list of read IDs""" 
	
	all_reads = all_reads.set_index('fragment_id') 
	
	#create df with only the entries in read_list
	new_df = all_reads.loc[read_list]

	#reset index on the new dataframe 
	new_df['fragment_id'] = new_df.index 
	new_df = new_df.reset_index(drop=True)
	del new_df['Unnamed: 0']

	return new_df

def pipelineDataFrame(read_list, all_reads): 
	"""Creates a dataframe of a subset of  reads using a list of read IDs""" 
	
	all_reads = all_reads.set_index('ReadID') 
	
	#create df with only the entries in read_list
	new_df = all_reads.loc[read_list]

	#reset index on the new dataframe 
	new_df['ReadID'] = new_df.index 
	new_df = new_df.reset_index(drop=True)

	return new_df


def findMissed(all_reads, detected_Vreads): 

	"""Looks up the integrations missed by the pipeline enabling subsequent determination of what kind of reads we missed""" 

	#find the maximum number of integrations that could be found from the reads if the pipeline was perfect 
	detectable_hPos = hPosList(all_reads) 
	print('\nMax number of integrations detectable from reads: '+str(len(detectable_hPos)), flush = True) 

	#look at how many of these integrations we cover with our results from the pipeline 
	#make a dataframe of the correct reads 
	correct_df = readDataFrame(detected_Vreads, all_reads) 
	#find the number of integrations that the pipeline detected 
	pipeline_hPos = hPosList(correct_df) 
	print("Number of integrations detected by the pipeline: "+str(len(pipeline_hPos)),flush = True) 

	acc = (len(pipeline_hPos)/len(detectable_hPos))*100
	print("Of the integrations which were present in the reads, {:.2f}% could be detected".format(acc),flush = True) 
   
	#create a column of the missed integrations so we can evaluate them
	missed_hPos =  []
	for hPos in detectable_hPos: 
		if hPos not in pipeline_hPos: 
			missed_hPos.append(hPos) 


	return missed_hPos

def hPosList(df): 
	"""Function to format the hPos in the output files into something which we can work with""" 
	#list to store the hPos 
	hPos_list = []

	for i in range(len(df)): 
		values = str(df['hPos'][i][1:-1])
		values = values.split(", ")
		values = list(map(int,values))
		for value in values: 
			if value not in hPos_list:
				hPos_list.append(value)  
				
	return hPos_list 
	 
def evaluateMissed(missed_hPos,ints):
	"""Attempt to find out why the integrations were completely missed by the pipeline""" 
	#create a dictionary of properities of each of the integrations
	hPos_dict = ints.set_index('hPos').to_dict()
	hPos_keys = hPos_dict.get('rearrangement').keys()

	#check out the different properities to see if we can find any particular issues 
	miss_rearrange = [hPos_dict.get('rearrangement').get(missed) for missed in missed_hPos]
	all_rearrange = [hPos_dict.get('rearrangement').get(key) for key in hPos_keys]
	miss_deletion = [hPos_dict.get('deletion').get(missed) for missed in missed_hPos]
	all_deletion = [hPos_dict.get('deletion').get(key) for key in hPos_keys]
	miss_num_fragments = [hPos_dict.get('num_fragments').get(missed) for missed in missed_hPos]
	all_num_fragments = [hPos_dict.get('num_fragments').get(key) for key in hPos_keys]

	#get statistics on the different properties 
	#get statisitics on rearangements 
	print("\n{:.2f}% of the total integrations were rearranged".format(percentactual(all_rearrange)))
	print("{:.2f}% of the integrations not detected by the pipeline were rearraged".format(percentactual(miss_rearrange)))
	#get statistics on deletions 
	print("\n{:.2f}% of the total integrations had deletions".format(percentactual(all_deletion)))
	print("{:.2f}% of the integrations not detected by the pipeline had deletions\n".format(percentactual(miss_deletion)))
	#get the statistics on the number of fragments use to separate the fragments 
	print("For all of the integrations...") 
	assessFragments(all_num_fragments)
	print("\nFor integrations missed by the pipeline...") 
	assessFragments(miss_num_fragments)   


def assessMissedLength(missed_hPos, ints): 
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
	
	#save distribution of all integration lengths   
	"""plt.pyplot.figure(0)
	weights_all = np.ones_like(all_lengths)/float(len(all_lengths)) 	
	plt.pyplot.hist(all_lengths, weights = weights_all) 
	plt.pyplot.title("Histogram of the length of all integrations") 
	plt.pyplot.xlabel('length of integration')
	plt.pyplot.ylabel('Density')  
	plt.pyplot.savefig("evaluate_pipeline_output/integration_lengths.pdf")
	print("Distribution of all integrations saved to integration_lengths.pdf") 
"""
	#save distribution of the missed integration lengths
	missed_ints = [hPos_length.get(missed_int) for missed_int in missed_hPos]
	weights_missed = np.ones_like(missed_ints)/float(len(missed_ints)) 
	plt.pyplot.figure(1)  
	plt.pyplot.hist(missed_ints, weights = weights_missed)
	plt.pyplot.title("Histogram of the length of missed integrations") 
	plt.pyplot.xlabel('length of integration')
	plt.pyplot.ylabel('Density')  
	plt.pyplot.savefig("evaluate_pipeline_output/missed_integration_lengths.pdf")
	print("Distribution of all integrations saved to evaluate_pipeline_output/missed_integration_lengths.pdf")
	
def percentactual(values): 
	"""counts the number of actual variables in a list of values""" 

	#find the number of missed integrations which have been rearranged 
	count  = values.count(actual) 

	per_actual = (count/len(values))*100 

	return per_actual

def assessFragments(num_fragments): 
	"""Tells us what number of fragments the missed integrations were divided into""" 
	min_frag = min(num_fragments) 
	max_frag = max(num_fragments)


	for i in range(min_frag, max_frag+1):
		count = num_fragments.count(i)
		frag_freq = (count/len(num_fragments))*100 
		print("{:.2f}% of integrations were broken into ".format(frag_freq)+str(i)+" fragments")


def getOverlap(all_reads):
	"""Tells us what type of junctions the filter/filtered reads had (from the actual file not the pipeline)
Inormation corresponds to read ID""" 
	#get the requried columns 
	left_read= all_reads['left_read'].values
	right_read = all_reads['right_read'].values
	left_junc = all_reads['left_junc'].values
	right_junc = all_reads['right_junc'].values
	int_ID = all_reads['fragment_id'].values

	#column which says the type of overlaps 
	overlap_type = []
	

	for i in range(len(all_reads)):
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

	
	
def assessLength(all_reads, undetected_Vreads): #TODO this is questionable remove 
	"""Compare the number of viral base pairs in the missed reads and in alfg
l reads"""
	#get the required columns
	int_ID = all_reads['fragment_id'].values
	first_len = all_reads['left_read_amount'].values
	second_len = all_reads['right_read_amount'].values

	#find the length of virus in each read - excluding reads consisting of whole virus and no virus
	viral_len = [max(first_len[i],second_len[i]) for i in range(len(all_reads))]	
	
	#create a dictionay
	len_dict = dict(zip(int_ID, viral_len)) 
	#get the lengths for missed reads 
	missed_len = [len_dict.get(missed) for missed in undetected_Vreads] 
	
	#save distribution of length of the missed data 
	ax = sns.distplot(missed_len, hist=False)
	ax.set_title("Distribton of the length of missed integrations") 
	ax.set(xlabel='length of missed integration', ylabel = 'Density')  
	fig = ax.get_figure()
	fig.savefig("evaluate_pipeline_output/missedread_len_distplot.pdf")
	print("Distribution of the length of the missed reads saved to evaluate_pipeline_output/evaluate_pipeline_output/missedread_len_distplot.pdf")


	
#def assessFalseReads #TODO look at the integrations which were falsely detected by the pipeline 
#look at the junctions in them 
#were they close to an actual integration event? 

def assessOverlap(all_reads, undetected_Vreads): 
	"""Looks at what type of overlap the reads missed by the pipeline have. Takes a dataframe of all of the reads and a list of the IDs of the missed reads""" 
	#dictionary of overlap types for each read 
	overlap_dict = getOverlap(all_reads)
	
	#column of the overlap of the undetected_Vreads
	overlap = [overlap_dict.get(missed) for missed in undetected_Vreads] 
	
	#find the instances of the gap types 
	gap = overlap.count('gap')
	overlap_ = overlap.count('overlap')
	none = overlap.count('none') #clean 
	other = overlap.count('nan') 

	#print the instances of each 
	gapFreq = (gap/len(undetected_Vreads))*100 
	print(str(gapFreq)+"\n% of missed reads had a gap") 
	overlapFreq = (overlap_/len(undetected_Vreads))*100
	print(str(overlapFreq)+"% of missed reads had an overlap") 
	cleanFreq = (none/len(undetected_Vreads))*100
	print(str(cleanFreq)+"% of missed reads had clean junctions") 
	otherFreq =(other/len(undetected_Vreads))*100
	print(str(otherFreq)+"%of missed reads had a junction other than those above") 


def pipelineMissed(missed_df,all_reads): 
	"""looks at which integrations cuase false negative reads"""

	#create a list of all possible hPos 
	all_hPos = hPosList(all_reads) 

	#create a list of the hpos in the missed viral reads 
	missr_hPos = []
	for i in range(len(missed_df)): 
		values = str(missed_df['hPos'][i][1:-1])
		values = values.split(", ")
		values = list(map(int,values))
		missr_hPos.append(values[0]) 

	#iterate through all_hPos counting the number of observations of each
	count_list = pd.Index(missr_hPos).value_counts(dropna = False) 
	count_list = pd.DataFrame({"hPos":count_list.index,"Count":count_list})
	
	#make a bar graph 
	fig, ax = plt.pyplot.subplots()  
	ax = sns.barplot(x = "hPos",y="Count", data = count_list)
	ax.set_ylabel("Frequency") 
	ax.tick_params(axis='both',which ='major', labelsize =3)
	ax.set_title("Frequency of each integration amongst false negative reads") 
	fig = ax.get_figure()
	fig.savefig("evaluate_pipeline_output/missedread_hPosfrequency.pdf")
	print("hPos frequencies saved to missed_read_hPosfrequency.pdf", flush=True) 

def missedLocation(missed_df): 
	"""plot the location of the location of the false positive reads in the genome - other homologous parts we dont know about?""" 
	
	#create a list of the read IDs 
	missed_IDs = missed_df["fragment_id"]

	#create a list of the coordinates spanning the read corresponding to missed_IDs
	missed_x = []
	missed_y = []

	for i in range(len(missed_df)): 
		if missed_df["left_read"][i] == "chimeric":
			values = str(missed_df["left_read_coor"][i][1:-1])
		else: 
			values = str(missed_df["right_read_coor"][i][1:-1])
		values = values.split(", ")
		values = list(map(int,values))
		missed_x.append(values[0])
		missed_y.append(values[1]) 

	#visualise as a boxplot 
	fig, ax = plt.pyplot.subplots()
	data = pd.DataFrame({"x":missed_x,"y":missed_y})
	ax = sns.scatterplot(x="x", y="y", data=data) 
	ax.set_ylabel("Reads") 
	ax.set_label("Location")
	ax.set_title("Location of false positive reads in the genome")
	ax.tick_params(axis="both",which="major", labelsize = 7) 
	fig = ax.get_figure()
	fig.savefig("evaluate_pipeline_output/missedread_location.pdf")

def readLength(missed_df): 
	"""Look at the length of viral DNA in false negative reads""" 
	
	#get the amount of viral DNA in the chimeric read 
	missed_viral = [] 
	
	for i in range(len(missed_df)): 
		if missed_df["left_read"][i] == "chimeric": 
			missed_viral.append(missed_df["left_read_amount"][i])
		else: 
			missed_viral.append(missed_df["right_read_amount"][i])

	#visualise as a distribution 
	fig, ax = plt.pyplot.subplots()
	ax = sns.distplot(missed_viral) #TODO normalise? 
	ax.set_xlabel("Viral DNA (bp)")
	ax.set_ylabel("Density") 
	ax.set_title("Distribution of the amount of viral DNA in false negative reads")
	fg = ax.get_figure()
	fig.savefig("evaluate_pipeline_output/missedread_viralDNA.pdf")
	print("Amount of viral DNA in reads missed by the pipeline saved to evaluate_pipeline_output/missedread_viralDNA.pdf",flush = True) 
	
	

def compareOverlap(all_reads, pipe_ints): #TODO finish this 
	"""Compares whether the pipeline has correctly predicted the gap type of the read. Takes list of predicted"""

	#get the IDs of the pipeline reads
	pipe_reads = pipe_ints['ReadID'].values 

	#get the overlap type of the pipeline reads
	pipe_overlap = pipe_ints['Type'].values
 
	#create a dictionary of the overlap types for each read 
	overlap_dict = getOverlap(all_reads)
	
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
