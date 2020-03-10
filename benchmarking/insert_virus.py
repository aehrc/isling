#### simulates viral insertions ####

#### written by Suzanne Scott (suzanne.scott@csiro.au) and Susie Grigson (susie.grigson@flinders.edu) ####

# types of possible insertions:
## whole viral genome (possibility for reverse orientation)
## portion of viral genome (possibility for reverse orientation)
## n sequential portions of viral genome, rearranged (with possibility that some are reversed)
## n non-sequential (with a gap in between) portions of viral genome, rearranged (with possibility that some are reversed) (ie deletion)
## portion of viral genome divided into n sequential portions of viral genome, rearranged 
## poriton of viral genome divided into n non-sequential portions of viral genome and rearranged 

#at each insertion type, overlaps possible are gap, overlap, none
## if gap, insert small number of random bases between host and virus
## if overlap, take small (<=10) number of bases and look for homology, do insertion there
## if none, there is a 'clean' junction between virus and host

# reports parameters of integrations:
	#location in host genome (chr, start, stop)
	#part of viral genome inserted (virus, start, stop)
	#integration type (whole, portion, rearrangement, etc)
	#overlaps/gaps at each junction 
# reports all integration sites relative to the host genome, independent of other integrations
# intergrations are stored internally though the Integration class

###import libraries
from Bio import SeqIO
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
import pandas as pd
import argparse
import sys
import os
import numpy as np
print("STARTING", flush=True) 
###
max_attempts = 50 #maximum number of times to try to place an integration site 

### main
def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='simulate viral integrations')
	parser.add_argument('--host', help='host fasta file', required = True, type=str)
	parser.add_argument('--virus', help = 'virus fasta file', required = True, type=str)
	parser.add_argument('--ints', help = 'output fasta file', required = True, type=str)
	parser.add_argument('--locs', help = 'output text with viral integrations', required = True, type=str)
	parser.add_argument('--ints_host', help = 'output csv with integration locations in host genome', required = True, type=str)
	parser.add_argument('--int_num', help = 'number of integrations to be carried out', required=True, type=int)
	parser.add_argument('--fasta', help = 'output fasta of integrated host genome', required = False, type=str)
	parser.add_argument('--sep', help = 'integrations must be separated by this many bases [20]', required=False, default=20, type=int)
	parser.add_argument('--min_len', help = 'minimum length of integerations [50]', required=False, default=50, type=int)
	parser.add_argument('--set_len', help = 'use to get integrations of a specific length', required=False, default=0, type=int)
	parser.add_argument('--epi_num', help = 'number of episomes [0]', required=False, default=0, type=int)
	parser.add_argument('--int_portion', help = 'specify is a particular type of integrations are wanted: whole, portion, both [both]', required=False, default='both', type=str)
	parser.add_argument('--int_deletion', help = 'specify if deletions should be conducted on none, some or all integrations [some]', required=False, default='rand', type=str)
	parser.add_argument('--int_rearrange', help = 'specify if rearrangements should be conducted on none, some or all integrations [some]', required=False, default='rand', type=str)
	parser.add_argument('--set_junc', help = 'specify is a particular type of junctions are wanted from: clean, gap, overlap or rand [rand]', required=False, default='rand', type=str)
	parser.add_argument('--seed', help = 'seed for random number generator', required=False, default=1, type=int)
	parser.add_argument('--verbose', help = 'display extra output for debugging', required=False, default=False, type=bool)
	args = parser.parse_args()
	
	#read host fasta - use index which doesn't load sequences into memory because host is large genome
	if checkFastaExists(args.host):
		host = SeqIO.index(args.host, 'fasta', alphabet=unambiguous_dna)
	else:
		raise OSError("Could not open host fasta")
	
	#read virus fasta -  make dictionary in memory of sequences
	if checkFastaExists(args.virus):
		virus = SeqIO.to_dict(SeqIO.parse(args.virus, 'fasta', alphabet=unambiguous_dna))
	else:
		raise OSError("Could not open virus fasta")

	# check that input arguments make sense
	# check that minimum integration size is longer than all the virus sequences
	if (args.set_len == 0) and not(all([args.min_len < len(seq) for seq in virus.values()])):
		raise ValueError(f"Minimum integration length (--min_len) must be longer than all virus sequences. \
				   viral sequence lengths: {[len(seq) for seq in virus.values()]}")

	# check that set integration size is longer than all the virus sequences
	if (args.set_len != 0) and (not(all([args.set_len < len(seq) for seq in virus.values()]))):
		raise ValueError(f"Minimum integration length (--set_len) must be longer than all virus sequences. \
				   viral sequence lengths: {[len(seq) for seq in virus.values()]}")

	# check that minimum separation allows for the desired number of integrations
	# assume that number of bases that are ruled out for each integration is 2 x sep
	if np.sum([len(seq) for seq in host.values()]) < (args.sep * 2):
		raise ValueError("Integrations are likely to clash. Specify fewer integrations or a shorter separation")

	# check the specified type of junction is valid
	if args.set_junc not in ["clean", "gap", "overlap", "rand"]:
		raise ValueError("junction type (--set_junc) must be one of \"clean\", \"gap\" or \"overlap\"")

	# check the specified type of integration is valid
	if args.int_portion not in ["whole", "portion", "both"]:
		raise ValueError("Not a valid type of integration: \"--int_portion\" should be one of \"whole\", \"portion\", \"both\"")
	if args.int_deletion not in ["all", "some", "none", "rand"]:
		raise ValueError("not a valid type of integration: \"--int_deletion\" should be one of \"all\", \"some\", \"none\"")
	if args.int_rearrange not in ["all", "some", "none", "rand"]:
		raise ValueError("not a valid type of integration: \"--int_rearrange\" should be one of \"all\", \"some\", \"none\"")

	# get integration types from combination of int_portion, int_deletion and int_rearrange
	int_types = parseIntTypes(args)

	
	# if the user has specified a set length, set the minimum length to be the same
	if args.set_len != 0:
		args.min_len = args.set_len

	#set random seed
	np.random.seed(args.seed)

	# get length to input into integration functions:
	# a string with either 'min' or 'set' and the corresponding value, separated by underscore
	if args.set_len == 0:
		int_len = "_".join(["min", str(args.min_len)])
	else:
		int_len = "_".join(["set", str(args.set_len)])

	#types of episomes 
	episome_types = [Episome.insertWhole, Episome.insertPortion, Episome.insertWholeRearrange, Episome.insertPortionRearrange, Episome.insertWholeDeletion, Episome.insertPortionDeletion] 
	
	#list to store integration objects
	host_ints = []
	#dictionary to  store sequences with integrations
	host_fasta = {key:value for key, value in zip(host.keys(), host.values())}
	
	#open output file for writing
	handle = open(args.locs, "w+")
	header = "\t".join(["hChr",
						"hPos",
						"inserted",
						"num_fragments",
						"virus",
						"oris",
						"rearrangement",
						"deletion",
						"vStart",
						"vStop",
						"leftJunction",
						"rightJunction",
						"breakpoints", 
						"vBases\n"
							])
	handle.write(header)
	
	#### PERFORM INTEGRATIONS ##### 
	
	#intialise after how many intergrations the number of integrations performed is reported to the user 
	int_report = 50
	
	print("\nNUMBER OF INTEGRATIONS TO DO: "+str(args.int_num),flush=True)
	
	#integration loop 
	counter = 0 # count number of iterations 

	while len(host_ints) < args.int_num: 

		#if a set length of integrations has been specified: 
		if args.set_len != 0: 

			host_ints, host_fasta = insertSetLength(host_fasta, virus, host_ints, handle, args.min_len, args.sep, args.set_len, args.set_junc)
		else: 
			#do random integrations if a specific type of integrations is not selected
			if args.int_type == 'rand': 
				int_type =  np.random.choice(insertion_types)

			#apply integration
			if args.verbose is True:	
				print(f"trying integration type: {int_type}, integrations already done: {[repr(int) for int in host_ints]}")
			host_ints, host_fasta = int_type(host_fasta, virus, host_ints, handle, args.min_len, args.sep, args.set_junc)
			

		#count the number of integrations applied 
		counter += 1  

		#periodically print how many integrations have been completed
		if counter % int_report == 0: 
			if args.verbose is True:
				print(str(counter) +" integrations attempted...", flush = True)

	print("NUMBER OF INTEGRATIONS DONE: "+str(len(host_ints)),flush=True)
	
	#episome loop 	
	print("\nNUMBER OF EPISOMES: "+str(args.epi_num))
	for i in range(0,args.epi_num): 
		rand_int = np.random.randint(0,2)	
		name = "episome "+str(i+1)
		host_fasta = episome_types[rand_int](virus, args.min_len, host_fasta, name)  
			
	print("\n***INTEGRATIONS COMPLETE***",flush=True)
	print(host_fasta,flush=True)
	handle.close()
	
	#save statistics on the integration 
	stats = Statistics.saveStats(host_ints)
	with open(args.ints_host, 'w') as handle:
		stats.to_csv(handle,sep='\t')
	
	#save integrated host sequence 
	with open(args.ints, 'w') as handle: 
    		SeqIO.write(host_fasta.values(), handle, 'fasta')
    		handle.close()
    		print("\nIntegrated host saved as "+str(args.ints),flush=True)
    		print("Details of sequences integrated saved as "+str(args.locs),flush=True)
    		print("Details of where intgrations lie in host sequence saved as "+str(args.ints_host),flush=True)

def insertWholeVirus(host, viruses, int_list, filehandle, min_len, sep, set_junc):
	"""Inserts whole viral genome into host genome"""

	#get positions of all current integrations
	currentPos = Statistics.integratedIndices(int_list)

	#make sure that new integration is not within args.sep bases of any current integrations
	attempts = 0
 
	while True:
		#get one viral chunk
		currentInt = Integration(host, set_junc)
		currentInt.addFragment(viruses, min_len, part = "whole")
		attempts += 1
		#check that all integrations are sep away from new integration
		if all([ abs(currentInt.hPos-i) > sep for i in currentPos ]):
			break
		#if we've made a lot of attempts, give up
		elif attempts > max_attempts:
			return int_list, host
	
	#do integration
	host, status = currentInt.doIntegration(host, int_list,filehandle)
	#only save if integration was successful 
	if status == True:
	
		#write to output fileh
		currentInt.writeIntegration(filehandle)
	
		#append to int_list
		int_list.append(currentInt)
	
	return int_list, host

def insertSetLength(host, viruses, int_list, filehandle, min_len, sep, set_len, set_junc):
	"""Inserts virus of a specified length. Not included in integration loop but useful for ATAY's work"""

	if set_len == 0: 
		raise ValueError("Specifiy a set integration size")
	
	#get positions of all current integrations
	currentPos = Statistics.integratedIndices(int_list)

	#make sure that new integration is not within args.sep bases of any current integrations
	attempts = 0
 
	while True:
		#get one viral chunk
		currentInt = Integration(host, set_junc)
		currentInt.addFragment(viruses, min_len, part = set_len)
		attempts += 1
		#check that all integrations are sep away from new integration
		if all([ abs(currentInt.hPos-i) > sep for i in currentPos ]):
			break
		#if we've made a lot of attempts, give up
		elif attempts > max_attempts:
			return int_list, host
	
	#do integration
	host, status = currentInt.doIntegration(host, int_list,filehandle)
	#only save if integration was successful 
	if status == True:
	
		#write to output filehandle
		currentInt.writeIntegration(filehandle)
	
		#append to int_list
		int_list.append(currentInt)
	
	return int_list, host

def insertViralPortion(host, viruses, int_list, filehandle, min_len, sep, set_junc):
	"""Inserts portion of viral DNA into host genome"""
	
	#get positions of all current integrations
	currentPos = Statistics.integratedIndices(int_list)
	
	#make sure that new integration is not within args.sep bases of any current integrations
	attempts = 0
	while True:
		#get one viral chunk 
		currentInt = Integration(host, set_junc)
		currentInt.addFragment(viruses, min_len, part = "rand")
		attempts += 1
		
		#check that all integrations are sep away from new integration
		if all([ abs(currentInt.hPos-i) > sep for i in currentPos ]):
			break
		#if we've made a lot of attempts, give up
		elif attempts > max_attempts:
			return int_list, host
			
	#do integration 
	host, status  = currentInt.doIntegration(host, int_list,filehandle)
	
	#only save if integration was successful 
	if status == True: 
		
		#write to output file
		currentInt.writeIntegration(filehandle)
	
		#append to int_list
		int_list.append(currentInt)
	
	return int_list, host

def insertWholeRearrange(host, viruses, int_list, filehandle, min_len,sep, set_junc):
	""" Inserts a single portion of viral DNA with n rearrangements """
	
	#get positions of all current integrations
	currentPos = Statistics.integratedIndices(int_list)
	
	#make sure that new integration is not within args.sep bases of any current integrations
	attempts = 0
	while True:
		#get one viral chunk
		currentInt = Integration(host, set_junc)
		currentInt.addRearrange(viruses, min_len, part = "whole")
		attempts += 1
		
		#check that all integrations are args.sep away from new integration
		if all([ abs(currentInt.hPos-i) > sep for i in currentPos ]):
			break
		#if we've made a lot of attempts, give up
		elif attempts > max_attempts:
			return int_list, host
			
	#do integration
	host, status  = currentInt.doIntegration(host, int_list,filehandle)
	
	#only save if integration was successful 
	if status == True: 
	
		#write to output file
		currentInt.writeIntegration(filehandle)
	
		#append to int_list
		int_list.append(currentInt)

	return int_list, host

def insertWithDeletion(host, viruses, int_list, filehandle, min_len,sep, set_junc):
	""" Inserts n portions of viral DNA into host genome"""
	#get positions of all current integrations 
	currentPos = Statistics.integratedIndices(int_list)

	#make sure that new integration is not within args.sep bases of any current integrations
	attempts = 0
	while True:
		#get one viral chunk
		currentInt = Integration(host, set_junc)
		currentInt.addDeletion(viruses, min_len,part = "whole")
		attempts += 1
		#check that all integrations are args.sep away from new integration
		if all([ abs(currentInt.hPos-i) > sep for i in currentPos ]):
			break
		#if we've made a lot of attempts, give up
		elif attempts > max_attempts:
			return int_list, host
	
	#do integration
	host, status = currentInt.doIntegration(host,int_list,filehandle)
	
	#only save if integration was successful 
	if status == True:
		
		#write to output file
		currentInt.writeIntegration(filehandle)
	
		#append to int_list
		int_list.append(currentInt)
	
	return int_list, host

def insertPortionRearrange(host, viruses, int_list, filehandle, min_len,sep, set_junc): 
	"""Rearranges a portion of virus""" 

	#get positions of all current integrations
	currentPos = Statistics.integratedIndices(int_list)
	
	#make sure that new integration is not within args.sep bases of any current integrations
	attempts = 0
	while True:
		#get one viral chunk
		currentInt = Integration(host, set_junc)
		currentInt.addRearrange(viruses, min_len, part = "rand")
		attempts += 1
		print(currentInt)
		
		#check that all integrations are args.sep away from new integration
		if all([ abs(currentInt.hPos-i) > sep for i in currentPos ]):
			break
		#if we've made a lot of attempts, give up
		elif attempts > max_attempts:
			return int_list, host
			
	#do integration
	host, status  = currentInt.doIntegration(host, int_list,filehandle)
	
	#only save if integration was successful 
	if status == True: 
	
		#write to output file
		currentInt.writeIntegration(filehandle)

	
		#append to int_list
		int_list.append(currentInt)
	
	return int_list, host

def insertPortionDeletion(host, viruses, int_list, filehandle, min_len,sep, set_junc):
	"""Takes a segment of the virus and then deletes a segment from it and inserts the resulting portion. This have the potentional to be short - useful for simulating short reads later""" 

	#get positions of all current integrations 
	currentPos = Statistics.integratedIndices(int_list)

	#make sure that new integration is not within args.sep bases of any current integrations
	attempts = 0
	while True:
		#get one viral chunk
		currentInt = Integration(host, set_junc)
		currentInt.addDeletion(viruses, min_len, part = "rand")
		attempts += 1
		#check that all integrations are args.sep away from new integration
		if all([ abs(currentInt.hPos-i) > sep for i in currentPos ]):
			break
		#if we've made a lot of attempts, give up
		elif attempts > max_attempts:
			return int_list, host
	
	
	#do integration
	host, status = currentInt.doIntegration(host,int_list,filehandle)
	
	#only save if integration was successful 
	if status == True:
		
		#write to output file
		currentInt.writeIntegration(filehandle)
	
		#append to int_list
		int_list.append(currentInt)
	
	
	return int_list, host


def checkFastaExists(file):
	#check file exists
	exists = os.path.isfile(file)
	if not(exists):
		return False
	#check extension
	prefix = file.split(".")[-1]
	if prefix:
		if not((prefix == "fa") or (prefix == "fasta") or (prefix == "fna")):
			return False
	return True
	
def parseIntTypes(args):
		# for convenience, if we want integrations of only one type to be "all", should be able to specify just that one type
	# but this creates conflicts between "all" and "some"
	# resolve this by using "rand" as default and changing values depending on inputs
	
	# if one all, assume the other should be none
	if (args.int_deletion == "all") and (args.int_rearrange == "rand"):
		args.int_rearrange = "none"
	if (args.int_rearrange == "all") and (args.int_deletion == "rand"):
		args.int_deletion = "none"
	# if one is some, assume the other should be some
	if (args.int_deletion == "some") and (args.int_rearrange == "rand"):
		args.int_rearrange = "some"
	if (args.int_rearrange == "some") and (args.int_deletion == "rand"):
		args.int_deletion = "some"
	# if one is none, the other should be some
	if (args.int_deletion == "none") and (args.int_rearrange == "rand"):
		args.int_rearrange = "some"
	if (args.int_rearrange == "none") and (args.int_deletion == "rand"):
		args.int_deletion = "some"
	# if neither specified, both should be some
	if (args.int_rearrange == "rand") and (args.int_deletion == "rand"):
		args.int_rearrange = "some"
		args.int_deletion = "some"

	#integration types are: [insertWholeVirus, insertViralPortion, insertWholeRearrange, insertWithDeletion, insertPortionRearrange, insertPortionDeletion]
	
	# only whole integrations
	if args.int_portion == "whole": 
	
		if (args.int_deletion == "all") and (args.int_rearrange == "none"):
			int_types = [insertWithDeletion]
		elif (args.int_deletion == "none") and (args.int_rearrange == "all"):
			int_types = [insertWholeRearrange]
		elif (args.int_deletion == "some") and (args.int_rearrange == "some"):
			int_types = [insertWholeVirus, insertWholeRearrange, insertWithDeletion]
		elif (args.int_deletion == "some") and (args.int_rearrange == "none"):
			int_types = [insertWholeVirus, insertWithDeletion]
		elif (args.int_deletion == "none") and (args.int_rearrange == "some"):
			int_types = [insertWholeVirus, insertWholeRearrange]
		elif (args.int_deletion == "none") and (args.int_rearrange == "none"):
			int_types = [insertWholeVirus]
		else:
			raise ValueError("combination of int_deletion and int_rearrange is invalid")
			
	elif args.int_portion == "portion": #only partial integrations	
		if (args.int_deletion == "all") and (args.int_rearrange == "none"):
			int_types = [insertPortionDeletion]
		elif (args.int_deletion == "none") and (args.int_rearrange == "all"):
			int_types = [insertPortionRearrange]
		elif (args.int_deletion == "some") and (args.int_rearrange == "some"):
			int_types = [insertViralPortion, insertPortionRearrange, insertPortionDeletion]
		elif (args.int_deletion == "some") and (args.int_rearrange == "none"):
			int_types = [insertViralPortion, insertPortionDeletion]
		elif (args.int_deletion == "none") and (args.int_rearrange == "some"):
			int_types = [insertViralPortion, insertPortionRearrange]
		elif (args.int_deletion == "none") and (args.int_rearrange == "none"):
			int_types = [insertViralPortion]
		else:
			raise ValueError("combination of int_deletion and int_rearrange is invalid")

	# both whole and partial integrations
	else:
		if (args.int_deletion == "all") and (args.int_rearrange == "none"):
			int_types = [insertWithDeletion, insertPortionDeletion]
		elif (args.int_deletion == "none") and (args.int_rearrange == "all"):
			int_types = [insertWholeRearrange,insertPortionRearrange]
		elif (args.int_deletion == "some") and (args.int_rearrange == "some"):
			int_types = [insertWholeVirus, insertViralPortion, insertWholeRearrange, insertWithDeletion, insertPortionRearrange, insertPortionDeletion]
		elif (args.int_deletion == "some") and (args.int_rearrange == "none"):
			int_types = [insertWholeVirus, insertViralPortion, insertWithDeletion, insertPortionDeletion]
		elif (args.int_deletion == "none") and (args.int_rearrange == "some"):
			int_types = [insertWholeVirus, insertViralPortion, insertWholeRearrange, insertPortionRearrange]
		elif (args.int_deletion == "none") and (args.int_rearrange == "none"):
			int_types = [insertWholeVirus, insertViralPortion]
		else:
			raise ValueError("combination of int_deletion and int_rearrange is invalid")
			
	if args.verbose is True:
		print(f"integration types will be: {int_types}")
			
	return int_types
	
class Integrations:
	"""
	
	"""
	
class Integration:
	""" 
	Class to store the properties of integrations 
	"""
	
	def __init__(self, host, set_junc):
		"""
		initalize integration at a random place in the host genome
		only one integration, but may be of a rearranged chunk
		"""
		self.host = host
		self.chr = np.random.choice(list(host.keys()))
		self.hPos = np.random.randint(1, len(host[self.chr].seq))
		self.fragments = 0
		
		#overlaps and gaps at each end of the integration
		#negative overlap means overlap (ie need to find where host and virus align)
		#positive overlap means gap (ie inserted bases that come from neither host nor virus)
		#zero overlap means clean junction
		#integration cannot not have negative overlap at both ends 


		junc_dict = {'clean': 0, 'gap': 1, 'overlap': 2}
		junction_types = ['clean', 'gap', 'overlap']

		if set_junc == 'rand': 
			junc_num = np.random.randint(0, len(junction_types)) 
		else: 
			if set_junc not in junction_types: 
				raise OSError("Not a valid type of junction")
			else: 
				junc_num = junc_dict.get(str(set_junc))
				
		this_junc = junction_types[junc_num] 

		#generate junctions 
		self.overlaps = self.addJunction(this_junc) 
		
		#record the type of junction for saving to file later
		self.junction = (self.convertJunction(self.overlaps[0]),
				self.convertJunction(self.overlaps[1]))

		#used when an integration with an overlap is integrated at a location other than hPos 
		self.newpoint = -1
		

	def addJunction(self, junc):
		""" sets the coordinates for a junction - clean, gap or overlap""" 
	

		#assign integration a clean junction
		if junc == 'clean':
			overlap1, overlap2 = (0,0)
		
		#assign integration a gap junction 
		elif junc == 'gap': 

			#generate a random number to decide whether there are gaps at one or both ends
			gap_num = np.random.randint(0,2) 

			#only have a gap at one end 
			if gap_num == 0: 
				
				#generate a random number to decide which end the gap goes on 
				gap_end = np.random.randint(0,2)
				
				#gap on the left 
				if gap_end == 0: 
					overlap1, overlap2 = (np.random.randint(1,10),0)

				#gap on the right 
				else: 
					overlap1, overlap2 = (0, np.random.randint(1,10))

			#have gaps at both ends 
			else: 
				overlap1, overlap2 = (np.random.randint(1,10), np.random.randint(1,10))	

		#assign integration an overlap juncion 
		elif junc == 'overlap':
		
			#generate a random number to decide which end the overlap ocurs at   
			end = np.random.randint(0,2)

			if end == 0: 
				overlap1, overlap2 = (np.random.randint(-10,0),0)
			else: 
				overlap1, overlap2 = (0,np.random.randint(-10,0)) 


		return overlap1, overlap2
 		 
	def addFragment(self, viruses, min_len, part = "rand"):
		"""
		add a viral fragment to this integration
		"""
		
		#check there isn't already a fragment
		if self.fragments > 0:
			print("Fragment has already been added!")
			return
			
		#add a simple fragment
		self.fragments = 1
		
		#get viral chunk
		self.chunk = ViralChunk(viruses, min_len, part)
	
		
	def addRearrange(self, viruses, min_len, part = "rand"):
		"""
		add a rearranged fragment
		"""
		
		#check there isn't already a fragment
		if self.fragments > 0:
			print("Fragment has already been added!")
			return
			
		#get number of fragments to rearrange into
		#draw from a poisson distribution with lambda = 1.5
		while True:
			n = int(np.random.poisson(1.5))
			if n > 1:
				break
		self.fragments = n
			
		#add a rearranged fragment
		self.chunk = ViralChunk(viruses, min_len, part)
		self.chunk.rearrange(n)

		#change the bases of the viral chunk to the rearranged fragment 
		new_chunk = ""
		for i in range(0,n): 
			portion = self.chunk.pieces.get(i).get('bases').seq
			new_chunk = new_chunk + portion
		self.chunk.bases = new_chunk
	
	def addDeletion(self, viruses, min_len,part = "rand"):
		"""
		add a fragment with a deletion in the middle
		always split into > 3 pieces and delete the middle
		"""
		#check there isn't already a fragment
		if self.fragments > 0:
			print("Fragment has already been added!")
			return
			
	#get number of fragments to rearrange into
		#draw from a poisson distribution with lambda = 2
		#make sure n > 3
		while True:
			n = int(np.random.poisson(2))
			if n > 3:
				break
		self.fragments = n
		
		#add a rearranged fragment
		self.chunk = ViralChunk(viruses, min_len, part)
		self.chunk.delete(n)
		self.fragments -= 1

		#create new chunk using the set of fragments
		new_chunk = ""

		#get the keys of the fragments after the deletion 
		frag_keys = self.chunk.pieces.keys()

		for i in frag_keys:
			portion = self.chunk.pieces.get(i).get('bases').seq
			new_chunk = new_chunk + portion 
			
		self.chunk.bases = new_chunk

	def createGap(self,size):
		"""
		Creates stretch of random DNA of length 'size' to insert as a gap 
		"""
		dna = ["a","g","c","t"]
		random_bases = ""
		for i in range(size):
			random_bases+=np.random.choice(dna)
		return random_bases

	def gapAdjust(self):
		"""
		Inserts random DNA where there is a gap junction 
		"""
		if self.overlaps[0]>0: #gap on left 
			gap = self.createGap(self.overlaps[0])
			self.chunk.bases = gap+self.chunk.bases
		if self.overlaps[1]>0: #gap on right 
			gap = self.createGap(self.overlaps[1])
			self.chunk.bases = self.chunk.bases+gap
		return
		
	def overlapPoint(self,left_site,right_site,filehandle):
		"""
		Finds point where there is homology between the human and viral sequences 
		"""
		overlap_point = 0
		#if no homologous region is found we skip the integration 
		if left_site ==-1 and right_site == -1: 
			self.chunk.bases = ""
			self.inserted = False 
			self.writeIntegration(filehandle)
		
		#if no homologous region found on left
		elif left_site ==-1: 
			overlap_point = right_site 
			
		#if no homologous region found on right 
		elif right_site ==-1: 
			overlap_point = left_site 
			
		#if both sides have homology 
		else: 
			if abs(left_site-self.hPos)<abs(right_site-self.hPos): 
				overlap_point = left_site
			else: 
				overlap_point = right_site

		#store the overlap point
		if self.inserted == True: 
			self.newpoint = overlap_point 


		return overlap_point  

	def createLeftOverlap(self,host,int_list,filehandle):
		"""
		Handles overlaps on the left of a viral chunk. Finds homolgous region closest to the integration site, making sure these homolougous regions are not from a previous integration
		"""
		
		#get the viral chunk 
		viral_chunk = self.chunk.bases
		#generate a list of locations where integrations have already occured 
		dont_integrate = self.dontIntegrate(int_list)


		#find sequence left of the integration 
		#convert the viral chunk to a seq object. Will already be a seq object if a deletion or rearrangement 
		if self.chunk.isRearranged == False and self.chunk.deletion == False:
			l_overlap = str(viral_chunk[:-self.overlaps[0]].seq)
		else: 
			l_overlap = str(viral_chunk[:-self.overlaps[0]])
			
		#give placeholder values to the homology sites left and right of the integration site 	
		left_site = -1
		right_site = -1

		#find homologous region on the left side 
		left_seq = host[self.chr][:self.hPos].seq
		while left_site <0: 
			left_search = str(left_seq).rfind(l_overlap)
			if left_search == -1:
				break 
			if left_search in dont_integrate: #check homology is not from a previous integration
				left_seq = left_seq[:left_search]  
			else:
				left_site= left_search 
				
		#find homologous region on the right side 
		right_seq = host[self.chr][self.hPos+self.overlaps[0]:].seq
		right_seqLen = len(right_seq)			
		while right_site >-1: #goes until no more sequence   
			right_search = str(right_seq).find(l_overlap)
			if right_search == -1:
				break 
			right_search += right_seqLen +(right_seqLen - len(right_seq))
			if right_search in dont_integrate: #check homology is not from a previous integration
				#start next search from after the homologous point from the integration was found 	
				right_seq = host[self.chr][right_search+self.overlaps[0]:]
			else:
				right_site = right_search
	
		overlap_point = self.overlapPoint(left_site,right_site,filehandle)
		int_start = overlap_point
		int_stop = overlap_point 


		if left_site != -1 and right_site != -1: 
		#only if an integration occurs does the overlapping region get removed 
			int_start =  0
			int_stop = 0
 
		return int_start,int_stop 
		
	def createRightOverlap(self,host,int_list,filehandle): 
		"""
		Handles overlaps on the right of a viral chunk. Finds homolgous region closest to the integration site, making sure these homolougous regions are not from a previous integration 
		""" 

		#get the viral chunk 
		viral_chunk = self.chunk.bases
		#generate a list of locations where integrations have already occured 
		dont_integrate = self.dontIntegrate(int_list)

		#find sequence right of the integration 
		#convert the viral chunk to a seq object. Will already be a seq object if a deletion or rearrangement 
		if self.chunk.isRearranged == False and self.chunk.deletion == False:
			r_overlap = str(viral_chunk[len(viral_chunk)+self.overlaps[1]:].seq)
		else: 
			r_overlap = str(viral_chunk[len(viral_chunk)+self.overlaps[1]:]) 

		#give placeholder values to the homology sites left and right of the integration site
		left_site = -1
		right_site = -1
	
		#find homologous region on the left side 
		left_seq = host[self.chr][:self.hPos].seq
		while left_site<0:
			left_search = left_seq.rfind(str(r_overlap))
			if left_search == -1:
				break 
			if left_search in dont_integrate: # check homology is not from a previous integration
				left_seq = left_seq[:left_search]  
			else: 
				left_site=left_search
	
		#find homologous region on the right side 
		right_seq = host[self.chr][self.hPos:].seq
		start_rseq = len(right_seq)
		right_spot = self.hPos+(start_rseq-len(right_seq))
			
		while right_site<0:
			right_search = str(right_seq).find(r_overlap) 
			if right_search == -1: 
				break 
			right_search += right_spot
			if right_search in dont_integrate: #check homology is not from a previous integration
				right_seq = host[self.chr][self.hPos:right_search]
			else: 
				right_site = right_search
				break

		overlap_point = self.overlapPoint(left_site,right_site,filehandle) 
		int_start = overlap_point
		int_stop = overlap_point
 
		if left_site != -1 and right_site != -1:
		#ensures only if there is homology the overlapping region is removed   
			int_stop = 0
			int_stop = 0
		return int_start,int_stop 
		
		
	def dontIntegrate(self,int_list): 
		"""
		Creates list of sites involved in viral integrations so that new integrations do not use these sites 
		"""

		dont_integrate = []
		
		stop_start = Statistics.adjustedStopStart(self,int_list)

		for i in range(len(stop_start)):
			c1, c2 = stop_start[i]

			for j in range(c1,c2): 
				dont_integrate.append(j)

			if int_list[i].overlaps[0] < 0:
				#add bases of left overlaps so any overlapping regions are not reused
				for k in range(0,abs(int_list[i].overlaps[0])+1): 
					dont_integrate.append(c1-k)
  
			if int_list[i].overlaps[1] < 0 : 
				#add bases of right overlaps so any overlapping regions are not reused 
				for k in range(0,abs(int_list[i].overlaps[1])+1): 
					dont_integrate.append(c2+k) 

		return dont_integrate	
		
			
	def doIntegration(self, host, int_list,filehandle):
		"""
		Inserts viral DNA (self.bases) at position self.hPos in host[self.chr]
		"""
		#start by assuming that integration will be successful 
		status = True
		
		#check there is already a fragment
		if self.fragments < 1:
			print("Add fragment before doing integration!")
			return
			
		#allows us to print whether or not the integration was successful 
		self.inserted = True 
		
		#need to account for previous integrations when inserting bases
		#get the number of bases previously added to this chromosome
		prevAdded = 0
		for int in int_list:
			if (int.chr == self.chr) & (int.hPos < self.hPos):
				prevAdded += int.numBases

		#keep track of how many bases added in this integration 
		self.numBases = self.overlaps[0] + self.overlaps[1] + len(self.chunk.bases)

		#use for inserting viral_chunk
		int_start = self.hPos 
		int_stop = self.hPos 

		#adjust seqences with gaps by adding random sequence
		self.gapAdjust()
		#list previous integrations 
		previousInt = Statistics.intList(self,int_list)

		#adjust left overlap 
		if self.overlaps[0]<0:
			(int_start,int_stop) = self.createLeftOverlap(host,int_list,filehandle)
			#remove homologous region so there is not double 
			self.chunk.bases = self.chunk.bases[-self.overlaps[0]:]
			if int_start != 0 and int_stop != 0: 
				int_start = int_start-self.overlaps[0]
			int_stop = int_start 
			 
		#adjust sequences with right overlap  		
		if self.overlaps[1]<0:
			(int_start,int_stop) = self.createRightOverlap(host,int_list,filehandle)
			int_stop = int_start 
			self.chunk.bases = self.chunk.bases[:len(self.chunk.bases)+self.overlaps[1]] 
	
			
		#If integration cannot be performed we stop 
		if int_stop == 0 and int_start == 0: 
			status = False
			
		host[self.chr] = host[self.chr][:int_start] + \
		 				"".join(self.chunk.bases) + \
		 				host[self.chr][int_stop:]
	
		
		#set the starting and stopping points of the performed integration 
		self.setStopStart(int_start,int_stop) 

		return host, status 


	def getOrisCoordsBases(self):
		"""
		return a lists of orientations, coordinates and bases for printing or writing to file
		"""
		if self.fragments > 0:
			pieces = self.chunk.pieces.keys()
			oris = [self.chunk.pieces[i]["ori"] for i in pieces]
			coords = [self.chunk.pieces[i]["coords"] for i in pieces]
			bases = [self.chunk.pieces[i]["bases"].seq for i in pieces]
	
			return oris, coords, bases
		else:
			return
	
	def writeIntegration(self, filehandle):
		#write information about the current integration to an output file
		#pass in an open filehandle for writing

		if self.fragments == 0:
			print("No fragments yet to write")
			return
			
		#get list of oris, breakpoints, bases
		(oris, coords, bases) = self.getOrisCoordsBases()
		#construct lines to write for current integration
		line = "\t".join([self.chr, 
							str(self.hPos),
							str(self.inserted), 
							str(self.fragments), 
							self.chunk.virus,
							"/".join(oris),
							str(self.chunk.isRearranged),
							str(self.chunk.deletion),
							str(self.chunk.start),
							str(self.chunk.stop),
							str(self.junction[0]),
							str(self.junction[1]),
							"/".join([f"{start}-{stop}" for (start, stop) in coords]),
							"/".join([str(x) for x in bases])+"\n"])
							
				
		filehandle.write(line)
			
		
	def convertJunction(self,value):
		"""
		Converts junction coordinates to names for saving to file
		""" 
		type = ""
		if value == 0: 
			type = "clean" 	
		elif value > 0: 
			type = "gap"
		elif value <0: 
			type = "overlap" 
		return type 
		
	def setStopStart(self,int_start,int_stop): 
		""" provides coordinates of where viral DNA starts and ends. hstart and hstop denote human genome"""

		if self.overlaps[0]>0: # ie there is a gap 
			self.hstart = int_start+self.overlaps[0]
		elif self.overlaps[0]<0: # ie there is an overlap
			self.hstart = int_start 
		else: 
			self.hstart = int_start 
			
		if self.overlaps[1]>0: # ie there is a gap
			self.hstop = int_stop-self.overlaps[1]
		elif self.overlaps[1]<0: #i e there is an overlap 
			self.hstop = int_stop
		else: 
			self.hstop = int_stop 
			
		self.hstop = self.hstop+len(self.chunk.bases) 

		
	def __repr__(self):
		if self.fragments == 0:
			return f"Integration on chromosome {self.chr}, position {self.hPos}"
		else:
			return f"Integration of {self.fragments} viral fragument on chromosome {self.chr}, position {self.hPos}.\n" \
			       f"Integrated virus: {repr(self.chunk)}"	
	def __str__(self):
			return f"Integration on chromosome {self.chr}, position {self.hPos}"		
		
class ViralChunk:
	"""
	Class with a chunk of virus
	Default is to get a random chunk (part = 'rand')
	use part = "whole" (or some other string) to get whole virus, or use part = int to get a chunk of lenth int
	"""
	
	#make viral chunk
	def __init__(self, viruses, min_len, part = "rand"):

		#get virus to integrate
		self.virus = np.random.choice(list(viruses.keys()))
		
		if min_len > len(viruses[self.virus].seq):
			raise ValueError("Viral genome is shorter than the minimum intgration size!")

		#if we want a random of chunk of a predetermined size (set_len) 
		if isinstance(part, int): 
			#get the length of the integration
			set_len = part 

			#can't use a set_len longer the virus
			if set_len > len(viruses[self.virus].seq): 
				raise ValueError("Set length for integrations is longer the viral sequence!")

			self.start = np.random.randint(0, len(viruses[self.virus].seq) - set_len)
			self.stop = self.start + set_len
		

		#if we want a random chunk of virus
		elif part == "rand":
			while True:
				self.start = np.random.randint(0, len(viruses[self.virus].seq) - min_len)
				self.stop = np.random.randint(self.start+1, len(viruses[self.virus].seq))
				if (self.stop - self.start) > min_len:
					break
		
		#if we want the whole virus
		else:
			self.start = 0
			self.stop = len(viruses[self.virus].seq)

		#give the integration an orientation 
		if np.random.uniform() > 0.5:
			self.ori = "f" #define orientation
			self.bases = viruses[self.virus][self.start:self.stop] #get bases to insert
		else:
			self.ori = "r" #define orientation
			self.bases = viruses[self.virus][self.start:self.stop].reverse_complement() #get bases to insert remove this one
		
		#construct dictionary with keys 'bases', 'ori' and 'coords'
		#use to keep track of order if 
		self.pieces = {0:{"bases":self.bases, "ori":self.ori, "coords":(self.start, self.stop)}}

		#store information about changes made to the chunk 
		self.isSplit = False #has chunk been split?
		self.isRearranged = False #has chunk been rearranged?
		self.deletion = False #does this chunk have a deletion? 
		
	def split(self, n):
		#split a part of a virus into n random parts
		#for later rearrangement or deletion

		if self.isSplit is True:
			print("Already split")
			return
		
		#need to check that we have enough bases to leave at least one base in each section
		if (self.stop - self.start) < n:
			print("Not enough bases to split")
			return
		
		#make orientations of each part the same as the original ori
		oris = [self.pieces[0]["ori"] for i in range(n)]
		
		#get n random coordinates to be breakpoints (plus original start and stop)
		while True:
			breaks = [self.start, self.stop] + list(np.random.randint(self.start+1, self.stop-1, n-1))
			# check that all coordinates are unique
			if len(set(breaks)) == len(breaks):
				break
				
		#sort breakpoints 
		breaks.sort()
				
		#get original bases of chunk to split
		if self.pieces[0]["ori"] == "f":
			bases_orig = self.pieces[0]["bases"]
		else:
			bases_orig = self.pieces[0]["bases"].reverse_complement()
		self.pieces = {} #clear dict to re-add pieces
		bases = []
		#split bases into pieces according to breakpoints
		for i in range(n):
			#get coords of this piece relative to bases we already extracted
			start = breaks[i] - breaks[0]
			stop = breaks[i+1] - breaks[0]
			#append bases  for this piece
			if oris[i] == 'f':
				bases.append(bases_orig[start:stop])
			else:
				bases.append(bases_orig[start:stop].reverse_complement())
		
		#convert oris, bases, breakpoints to dict
		self.pieces = { i:{"bases":bases[i], "ori":oris[i], "coords":(breaks[i], breaks[i+1])} for i in range(n)}
		
		# isSplit is now true
		self.isSplit = True
		
	def rearrange(self, n):
		#rearrange a chunk in n parts
		#pick n random breakpoints, and extract sequence of each
		#recombine in a random order, with equal probability of each
		#part being in forward or reverse order
		
		#check if chunk has already been split into parts
		if self.isSplit is False:
			self.split(n)
		#if has been split into a different number of parts, just use that number
		else:
			n = len(self.breakpoints)

		# check we were able to successfully split the chunk
		if self.isSplit is False:
			return

		#don't rearrange chunk if it's already been rearranged
		if self.isRearranged is True:
			raise ValueError("Chunk has already been rearranged!")
		
		#get new order, and make sure it's not the same as before
		order = np.arange(n)
		while True:
			np.random.shuffle(order)
			if not(all([i==j for i, j in zip(order, range(len(order)))])):
				break
		
		#replace keys with shuffled numbers
		self.pieces = {order[key]:value for key, value in self.pieces.items()}

		self.isRearranged = True #rearranged is now true
			
	def delete(self, n):
		#delete a middle from a viral chunk of n pieces
		
		#if there are less than 3 pieces, just use 3
		if n < 2:
			n = 3
	
		#check if chunk has already been split into parts
		if self.isSplit is False:
			self.split(n)

		#if has been split into a different number of parts, just use that number
		else:
			n = len(self.breakpoints)

		# check we were able to successfully split the chunk
		if self.isSplit is False:
			return

		#don't rearrange chunk if it's already been rearranged
		if self.isRearranged is True:
			raise ValueError("Chunk has already been rearranged!")

	
		#get piece to delete - use middle piece
		key_del = np.around((max(self.pieces.keys()) - min(self.pieces.keys()))/2)

		#replace keys with shuffled numbers
		del self.pieces[key_del]

		self.deletion = True #deletion is now true
		
	def __repr__(self):
		split = "split" if self.isSplit else "unsplit"
		rearranged = "rearranged" if self.isRearranged else "unrearranged"
		deletion = "deletion" if self.deletion else "no deletion"
		return f"{split}, {rearranged} viral chunk with {deletion} from virus {self.virus} with start and stop coordinates ({self.start}, {self.stop}), pieces [self.pieces[key] for key in self.pieces.keys()]."

class Statistics:
	"""
	Class of functions which give information on the integrations performed 
	
	"""
	
	def intList(self,int_list):
		"""Makes list of the sites of previously preformed integrations and adjusts for the bases added due to integration"""	
		previousInt = [int.hPos for int in int_list]
		for i in range(1,len(previousInt)):
			for j in range(0,i):
				if int_list[i].hPos < int_list[j].hPos:
					previousInt[j] = previousInt[j]+int_list[i].numBases		
		return previousInt
		
	def adjustedStopStart(self,int_list): 
		"""Makes a list of the coordinates of the viral integrations adjusted with insertions added"""
			
		intCoords = [(int.hstart,int.hstop) for int in int_list]

		for i in range(0,len(intCoords)):


			for j in range(0,i):

				#consider if an integration site was moved due to an overlap 
				if int_list[i].newpoint > -1: 
					i_site = intCoords[i][0]
				else: 
					i_site = intCoords[i][0] 

			
				#consider if an integration site was moved due to an overlap 
				if int_list[j].newpoint > -1: 
					j_site = intCoords[j][0] 
				else: 					
					j_site = intCoords[j][0] 
						

				if i_site < j_site:
					#if an integration has an overlap we don't want the indexing to be adjusted for these bases 
					if int_list[i].overlaps[0] < 0 or  int_list[i].overlaps[1] < 0:
						#handle overlap on left 

						if int_list[i].overlaps[0]<0: 
							shift = abs(int_list[i].overlaps[0])
						#handle overlap on the right 
						else: 
							shift = abs(int_list[i].overlaps[1]) 

						#start coordinate of the integration 
						new_coord1 = intCoords[j][0]+int_list[i].numBases 
						#stop coordinate of the integraton 
						new_coord2 = intCoords[j][1]+int_list[i].numBases 
						intCoords[j] = (new_coord1,new_coord2)
	 
						
					else: 
						#start coordinate of the integration 
						new_coord1 = intCoords[j][0]+int_list[i].numBases
						#stop coordinate of the integration  
						new_coord2 = intCoords[j][1]+int_list[i].numBases 
						intCoords[j] = (new_coord1,new_coord2)

		#adjust intCoords for differences in indexing
		for i in range(len(intCoords)): 
			(c1,c2) = intCoords[i]
			intCoords[i] = (c1, c2-1) 
 		
			
		return intCoords




	def integratedIndices(int_list): 
		"""Creates a lisit of sequence indexes which are of viral origin""" 
		#create list of viral coordinates 		
		viral_idx = []
		
		#number of integrations which have been performed
		num_ints = len(int_list)  
		#get start and stop coordinates 
		if num_ints != 0: 
			int_coords = Statistics.adjustedStopStart(int_list[1-num_ints],int_list) 
			for i in range(0,num_ints): 
				c1, c2 = int_coords[i]
				for j in range(c1,c2+1): 
					viral_idx.append(j) 
		return viral_idx 
		
		
	def saveStats(int_list):
		"""Save a csv file with the coordinates of the adjusted insertions so when we create artifical reads we can predict which reads contain viral DNA""" 
		
		
		num_ints = len(int_list) 
		stats_df = pd.DataFrame()
		
		#include hPos - is unique to each integration and allows subsequent identification of what occured 
		previousInt = [int.hPos for int in int_list]

		#include our junction types - used later to find what type of junction are present 
		left_overlap = [int.junction[0] for int in int_list]
		right_overlap = [int.junction[1] for int in int_list]

		#include the amount of junction - used later
		left_bases = [int.overlaps[0] for int in int_list]
		right_bases = [int.overlaps[1] for int in int_list]  
		#take absolute value of the number of base pairs in each list. Does not matter whether bases are homoloous or gaps - they aren't viral and shouldn't be included as part of the integration
		left_bases = [abs(bases) for bases in left_bases]
		right_bases = [abs(bases) for bases in right_bases] 
		
		#can only save stats if integrations have been performed 
		if num_ints != 0 :
			int_sites = Statistics.intList(int_list[1-num_ints],int_list)
			int_coords = Statistics.adjustedStopStart(int_list[1-num_ints],int_list) 	

			int_start = []
			int_stop = []
		
			for i in range(0,num_ints): 
				c1,c2 = int_coords[i]
				int_start.append(c1)
				int_stop.append(c2)
		
			stats_df = pd.DataFrame({"Integration sites":int_sites,"Start point":int_start,"Stop point": int_stop, "hPos":previousInt, "leftJunction":left_overlap, "rightJunction":right_overlap, "leftJunctionBases": left_bases, "rightJunctionBases": right_bases}) 
		else:
			print("NO SUCCESSFUL INTEGRATIONS WERE PERFORMED")
		
		return stats_df 

class Episome: 
	"""Class of functions used to create episomes to add to the resulting fasta file""" 

	
	def insertWhole(viruses, min_len, host, name):
		"""Adds a viral sequence to the output fasta file without modification"""
		#get a chunk of virus 
		chunk = ViralChunk(viruses, min_len, 'whole')
		
		#as the episome is circular it can be cut at any point 
		#randomly select the point at which the episome is 'cut' 
		cutP = np.random.randint(0,len(chunk.bases)) 

		#rejoin episome
		before_cut = chunk.bases[:cutP]
		after_cut = chunk.bases[cutP:]
		new_seq = after_cut + before_cut
		chunk.bases = new_seq
		
		#add episome to the host dictionary 
		entry = SeqRecord(Seq(str(chunk.bases.seq)), id = name, name = name, description = "") 
		host[name] = entry  
		
		return host 

	def insertPortion(viruses, min_len, host, name): 
		"""Adds a portion of the viral sequence to the output fasta file""" 
		#get a chunk of virus 
		chunk = ViralChunk(viruses, min_len, 'rand')

		#as the episome is circular it can be cut at any point 
		#randomly select the point at which the episome is 'cut' 
		cutP = np.random.randint(0,len(chunk.bases)) 

		#rejoin episome  
		before_cut = chunk.bases[:cutP]
		after_cut = chunk.bases[cutP:]
		new_seq = after_cut + before_cut
		chunk.bases = new_seq
		
		#add episome to the host dictionary 
		entry = SeqRecord(Seq(str(chunk.bases.seq)), id = name, name = name, description = "") 
		host[name] = entry  
		
		return host 

	def insertWholeRearrange(viruses, min_len, host, name): 
		"""Adds a rearranged viral sequece to the output fasta file""" 
		#need to split it then rearrange the split

		#get a chunk of virus 
		chunk = ViralChunk(viruses, min_len, 'whole')

		#as the episome is circular it can be cut at any point 
		#randomly select the point at which the episome is 'cut' 
		cutP = np.random.randint(0,len(chunk.bases)) 

		#get number of fragments to rearrange into episome
		#draw from a poisson distribution with lambda = 1.5
		while True:
			n = int(np.random.poisson(1.5))
			if n > 1:
				break

		chunk.rearrange(n) 

		#change the bases of the viral chunk to the rearranged fragment 
		new_chunk = ""
		for i in range(0,n): 
			portion = chunk.pieces.get(i).get('bases').seq
			new_chunk = new_chunk+portion
		chunk.bases = new_chunk

		#rejoin episome  
		before_cut = chunk.bases[:cutP]
		after_cut = chunk.bases[cutP:]
		new_seq = after_cut + before_cut
		chunk.bases = new_seq
		
		#add episome to the host dictionary 
		entry = SeqRecord(Seq(str(chunk.bases)), id = name, name = name, description = "") 
		host[name] = entry 
 
		return host

	def insertPortionRearrange(viruses, min_len, host, name):
		"""Adds a portion of a rearranged viral sequence to the output fasta file""" 

		#get a chunk of virus 
		chunk = ViralChunk(viruses, min_len, 'rand')

		#as the episome is circular it can be cut at any point 
		#randomly select the point at which the episome is 'cut' 
		cutP = np.random.randint(0,len(chunk.bases)) 

		#get number of fragments to rearrange into episome
		#draw from a poisson distribution with lambda = 1.5
		while True:
			n = int(np.random.poisson(1.5))
			if n > 1:
				break

		chunk.rearrange(n) 

		#change the bases of the viral chunk to the rearranged fragment 
		new_chunk = ""
		for i in range(0,n): 
			portion = chunk.pieces.get(i).get('bases').seq
			new_chunk = new_chunk+portion
		chunk.bases = new_chunk

		#rejoin episome  
		before_cut = chunk.bases[:cutP]
		after_cut = chunk.bases[cutP:]
		new_seq = after_cut + before_cut
		chunk.bases = new_seq
		
		#add episome to the host dictionary 
		entry = SeqRecord(Seq(str(chunk.bases)), id = name, name = name, description = "") 
		host[name] = entry 
 
		return host

	def insertWholeDeletion(viruses, min_len, host, name):
		"""Inserts a whole viral sequence with deletion to the output fasta file""" 
		#get an order for the fragments 
		 
		#get a chunk of virus 
		chunk = ViralChunk(viruses, min_len, 'whole')

		#randomly select the point at which the episome is 'cut' 
		cutP = np.random.randint(0,len(chunk.bases))

		#get number of fragments to rearrange into
		#draw from a poisson distribution with lambda = 2
		#make sure n > 3
		while True:
			n = int(np.random.poisson(2))
			if n > 3:
				break
		chunk.delete(n)
		new_chunk = ""

		#get the keys of the fragments after the deletion 
		frag_keys = chunk.pieces.keys()

		for i in frag_keys:
			portion = chunk.pieces.get(i).get('bases').seq
			new_chunk = new_chunk + portion 
			
		chunk.bases = new_chunk

		#rejoin episome  
		before_cut = chunk.bases[:cutP]
		after_cut = chunk.bases[cutP:]
		new_seq = after_cut + before_cut
		chunk.bases = new_seq
		
		#add episome to the host dictionary 
		entry = SeqRecord(Seq(str(chunk.bases)), id = name, name = name, description = "") 
		host[name] = entry 
 
		return host

	def insertPortionDeletion(viruses, min_len, host, name):
		"""Inserts a whole viral sequence with deletion to the output fasta file""" 
		#get an order for the fragments 
		 
		#get a chunk of virus 
		chunk = ViralChunk(viruses, min_len, 'rand')

		#randomly select the point at which the episome is 'cut' 
		cutP = np.random.randint(0,len(chunk.bases))

		#get number of fragments to rearrange into
		#draw from a poisson distribution with lambda = 2
		#make sure n > 3
		while True:
			n = int(np.random.poisson(2))
			if n > 3:
				break
		chunk.delete(n)
		new_chunk = ""

		#get the keys of the fragments after the deletion 
		frag_keys = chunk.pieces.keys()

		for i in frag_keys:
			portion = chunk.pieces.get(i).get('bases').seq
			new_chunk = new_chunk + portion 
			
		chunk.bases = new_chunk

		#rejoin episome  
		before_cut = chunk.bases[:cutP]
		after_cut = chunk.bases[cutP:]
		new_seq = after_cut + before_cut
		chunk.bases = new_seq
		
		#add episome to the host dictionary 
		entry = SeqRecord(Seq(str(chunk.bases)), id = name, name = name, description = "") 
		host[name] = entry 
 
		return host 
 

if __name__ == "__main__":
	main(sys.argv[1:])
