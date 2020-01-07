#### simulates viral insertions ####

# types of possible insertions:
## whole viral genome (possibility for reverse orientation)
## portion of viral genome (possibility for reverse orientation)
## n sequential portions of viral genome, rearranged (with possibility that some are reversed)
## n non-sequential (with a gap in between) portions of viral genome, rearranged (with possibility that some are reversed)

#at each insertion type, overlaps possible are gap, overlap, none
## if gap, insert small number of random bases between host and virus
## if overlap, take small (<5) number of bases and look for homology, do insertion there
## if none, there is a 'clean' junction between virus and host

# reports parameters of integrations:
	#location in host genome (chr, start, stop)
	#part of viral genome inserted (virus, start, stop)
	#integration type (whole, portion, rearrangement, etc)
	#overlaps/gaps at each junction - *yet to be implemented*
# reports all integration sites relative to the host genome, independent of other integrations
# so if there are two integrations on the same chromosome, both positions in host genome will be relative to reference

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
print("STARTING") 
###
max_attempts = 5 #maximum number of times to try to place an integration site

### main
def main(argv):
	#get arguments
	parser = argparse.ArgumentParser(description='simulate viral insertions')
	parser.add_argument('--host', help='host fasta file', required = True)
	parser.add_argument('--virus', help = 'virus fasta file', required = True)
	parser.add_argument('--ints', help = 'output fasta file', required = True)
	parser.add_argument('--locs', help = 'output text with viral integrations', required = True)
	parser.add_argument('--ints_host', help = 'output csv with integration locations in host genome', required = True)
	parser.add_argument('--int_num', help = 'number of integrations to be carried out', required=False, default=0)
	parser.add_argument('--fasta', help = 'output fasta of integrated host genome', required = False)
	parser.add_argument('--sep', help = 'integrations must be seperated by this many bases', required=False, default=5)
	parser.add_argument('--min_len', help = 'minimum length of integerations', required=False, default=50)
	parser.add_argument('--epi_num', help = 'number of episomes', required=False, default=0)
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
	
	#set random seed
	np.random.seed(27)

	#types of insertions
	insertion_types = [insertWholeVirus, insertViralPortion, insertWholeRearrange, insertWithDeletion, insertPortionRearrange, insertPortionDeletion]

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
						"vStart",
						"leftJunction",
						"rightJunction",
						"breakpoints", 
						"vBases\n"
							])
	handle.write(header)
	
	#### PERFORM INTEGRATION #####
	
	#intialise the required number of integrations 
	int_num = int(args.int_num)
	
	#intialise the minimum chunk size 
	min_len = int(args.min_len) 
	
	#intialise the separation from other integrations 
	sep = int(args.sep) 
	
	#intialise how many episomal sequences included in the outputted fasta file 
	epi_num = int(args.epi_num) 
	
	#intialise after how many intergrations the number of integrations performed is reported to the user 
	int_report = 50
	
	
	print("\nNUMBER OF INTEGRATIONS TO INSERT: "+str(int_num))
	print("All integrations performed as whole rearranges") #TODO update this to reflect changes as we go 
	
	#integration loop 
	for i in range(0,int_num):
		#rand_int =  np.random.randint(0,len(insertion_types))
		rand_int = 2
		host_ints, host_fasta = insertion_types[rand_int](host_fasta, virus, host_ints, handle, min_len, sep)  
		if i % int_report == 0 and i != 0: 
			print(str(i) +" integrations complete...")
	
	print("NUMBER OF INTEGRATIONS DONE: "+str(len(host_ints)))		
	print("\nNUMBER OF EPISOMES: "+str(epi_num))
	
	
	for i in range(0,epi_num): 
		rand_int = np.random.randint(0,2)	
		name = "episome "+str(i+1)
		host_fasta = episome_types[rand_int](virus, min_len, host_fasta, name)  
			
	print("\n***INTEGRATIONS COMPLETE***")
	print(host_fasta)
	handle.close()
	

	#save statistics on the integration 
	stats = Statistics.saveStats(host_ints)
	with open(args.ints_host, 'w') as handle:
		stats.to_csv(handle,sep='\t')
	
	#save integrated host sequence 
	with open(args.ints, 'w') as handle: 
    		SeqIO.write(host_fasta.values(), handle, 'fasta')
    		handle.close()
    		print("\nIntegrated host saved as "+str(args.ints))
    		print("Details of sequences integrated saved as "+str(args.locs))
    		print("Details of where intgrations lie in host sequence saved as "+str(args.ints_host))


def insertWholeVirus(host, viruses, int_list, filehandle, min_len, sep):
	"""Inserts whole viral genome into host genome"""
	
	#get positions of all current integrations
	currentPos = Statistics.integratedIndices(int_list)

	#make sure that new integration is not within args.sep bases of any current integrations
	attempts = 0

	while True:
		#get one viral chunk
		currentInt = Integration(host)
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
	
		#write to output file
		currentInt.writeIntegration(filehandle)
	
		#append to int_list
		int_list.append(currentInt)
	
	return int_list, host

def insertViralPortion(host, viruses, int_list, filehandle, min_len,sep):
	"""Inserts portion of viral DNA into host genome"""
	
	#get positions of all current integrations
	currentPos = Statistics.integratedIndices(int_list)
	
	
	#make sure that new integration is not within args.sep bases of any current integrations
	attempts = 0
	while True:
		#get one viral chunk
		
		currentInt = Integration(host)
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

def insertWholeRearrange(host, viruses, int_list, filehandle, min_len,sep):
	""" Inserts a single portion of viral DNA with n rearrangements """
	
	#get positions of all current integrations
	currentPos = Statistics.integratedIndices(int_list)
	
	#make sure that new integration is not within args.sep bases of any current integrations
	attempts = 0
	while True:
		#get one viral chunk
		currentInt = Integration(host)
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

def insertWithDeletion(host, viruses, int_list, filehandle, min_len,sep):
	""" Inserts n portions of viral DNA into host genome"""
	#get positions of all current integrations 
	currentPos = Statistics.integratedIndices(int_list)

	#make sure that new integration is not within args.sep bases of any current integrations
	attempts = 0
	while True:
		#get one viral chunk
		currentInt = Integration(host)
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

def insertPortionRearrange(host, viruses, int_list, filehandle, min_len,sep): 
	"""Rearranges a portion of virus""" 

	#get positions of all current integrations
	currentPos = Statistics.integratedIndices(int_list)
	
	#make sure that new integration is not within args.sep bases of any current integrations
	attempts = 0
	while True:
		#get one viral chunk
		currentInt = Integration(host)
		currentInt.addRearrange(viruses, min_len, part = "rand")
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

def insertPortionDeletion(host, viruses, int_list, filehandle, min_len,sep):
	"""Takes a segment of the virus and then deletes a segment from it and inserts the resulting portion. This have the potentional to be short - useful for simulating short reads later""" 

	#get positions of all current integrations 
	currentPos = Statistics.integratedIndices(int_list)

	#make sure that new integration is not within args.sep bases of any current integrations
	attempts = 0
	while True:
		#get one viral chunk
		currentInt = Integration(host)
		currentInt.addDeletion(viruses, min_len,part = "rand")
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
	
class Integrations:
	"""
	
	"""
	
class Integration:
	""" 
	Class to store the properties of integrations 
	"""
	
	def __init__(self, host):
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
		#integration cannot not have negative overlap at both ends - can be changed later** 

		while True: 
			self.overlaps = (np.random.randint(-10,10),np.random.randint(-10,10))
			if self.overlaps[0]<0 and self.overlaps[1]<0: 
				continue
			else: 
				break 
		
		#record the type of junction for saving to file later
		self.junction = (self.convertJunction(self.overlaps[0]),
				self.convertJunction(self.overlaps[1]))

		#used when an integration with an overlap is intgrated at a location other than hPos 
		self.newpoint = -1 
 		 
	def addFragment(self, viruses, min_len, part = "rand"):
		"""
		add a viral fragment to this integration
		"""
		
		#check there isn't alreafbdy a fragment
		if self.fragments > 0:
			print("Fragment has already been added!")
			return
			
		#add a simple fragment
		self.fragments = 1
		
		#get viral chunk
		self.chunk = ViralChunk(viruses, min_len, part)
	
	#def addFragmentRearrange()
		
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
			new_chunk = new_chunk+portion
		self.chunk.bases = new_chunk
	
	def addDeletion(self, viruses, min_len,part = "rand"):
		"""
		add a fragment with a deletion in the middle
		always split into >3 pieces and delete the middle
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
		Handles overlaps on the left of a viral chunk. Left and right are treated as different functions as different operations must be performed. 
		Works by finding closest regions of homology on left side of the randomly selected integration point. It checks if the homologous region is caused by an existing integration and concatenates the search range and attempts to find a homologous region again if the homology was caused by an existing integration. This is repeated for the right side of te integration point. The region closest to randomly selected integration point is then used to insert the viral chunk. 
		"""
 		
		viral_chunk = self.chunk.bases

		l_overlap = str(viral_chunk[:-self.overlaps[0]]) #sequence left of the integration site 
		left_site = -1
		right_site = -1
		dont_integrate = self.dontIntegrate(int_list)
				
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
				right_seq = host[self.chr][right_search:]
			else:
				right_site = right_search
		#print("right site: "+str(right_site)) #debugging remove 
		#print("left site: "+str(left_site)) #debugging remove 
		overlap_point = self.overlapPoint(left_site,right_site,filehandle)
		int_start = overlap_point
		int_stop = overlap_point 

		if left_site != -1 and right_site != -1: 
		#only if an integration occurs does the overlapping region get removed 
			int_start =  overlap_point + self.overlaps[0] 
		return int_start,int_stop 
		
	def createRightOverlap(self,host,int_list,filehandle): 
		"""
		Handles overlaps on the right of a viral chunk. Same as above with operations applicable to right end 
		""" 
		viral_chunk = self.chunk.bases
		dont_integrate = self.dontIntegrate(int_list)
		
		r_overlap = str(viral_chunk[len(viral_chunk)+self.overlaps[1]:]) #sequence right of the integration site 
		left_site = -1
		right_site = -1
	
		#find homologous region on the left side 
		left_seq = host[self.chr][:self.hPos].seq
		
		while left_site<0:
			left_search = str(left_seq).rfind(r_overlap)
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
			#right_spot = self.hPos+(start_rseq-len(right_seq)) #TODO remove 

		overlap_point = self.overlapPoint(left_site,right_site,filehandle) 

		int_start = overlap_point
		int_stop = overlap_point
		 
		if left_site != -1 and right_site != -1:
		#ensures only if there is homology the overlapping region is removed   
			int_stop = overlap_point-self.overlaps[1]
		return int_start,int_stop 
		
		
	def dontIntegrate(self,int_list): 
		"""
		Creates list of sites involved in viral integrations so that new integrations do not use these sites 
		"""
		#should this actually use the adjusted integration start and stops 
		dont_integrate = []
		
		stop_start = Statistics.adjustedStopStart(self,int_list)
		for i in range(len(stop_start)):
			c1, c2 = stop_start[i]
			for j in range(c1,c2+1): 
				dont_integrate.append(j) 
			
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
		#adjust sequences with overlap
		#list previous integrations 
		previousInt = Statistics.intList(self,int_list)

		#adjust left overlap 
		if self.overlaps[0]<0:
			(int_start,int_stop) = self.createLeftOverlap(host,int_list,filehandle)
			 
		#adjust sequences with right overlap  		
		if self.overlaps[1]<0:
			(int_start,int_stop) = self.createRightOverlap(host,int_list,filehandle)
			
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

		if self.overlaps[0]>0: #ie there is a gap 
			self.hstart = int_start+self.overlaps[0]
		elif self.overlaps[0]<0: #ie there is an overlap
			self.hstart = int_start-self.overlaps[0]
		else: 
			self.hstart = int_start 
			
		if self.overlaps[1]>0: #ie there is a gap
			self.hstop = int_stop-self.overlaps[1]
		elif self.overlaps[1]<0: #ie there is an overlap 
			self.hstop = int_stop+self.overlaps[1]
		else: 
			self.hstop = int_stop 
			
		self.hstop = self.hstop+len(self.chunk.bases) 

		
	def __str__(self):
		if self.fragments == 0:
			return f"Integration on chromosome {self.chr}, position {self.hPos}"
		elif self.fragments == 1:
			(oris, coords, bases) = getOrisCoordsBases(self)
			return f"Integration of {self.fragments} viral fragument on chromosome {self.chr}, position {self.hPos}.\n" \
				f"Integrated virus: {self.chunk.virus}, virus bases: {bases}, virus orientations: {oris}, coords {coords}"	
		else:
			(oris, coords, bases) = getOrisCoordsBases(self)
			return f"Integration of {self.fragments} viral fraguments on chromosome {self.chr}, position {self.hPos}.\n" \
				f"\tIntegrated virus: {self.chunk.virus},\n\tvirus bases: {bases}," \
				f"\n\tvirus orientations: {oris},\n\tvirus breakpoints: {coords}"			
		
class ViralChunk:
	"""
	Class with a chunk of virus
	Default is to get a random chunk (part = 'rand')
	use part = "whole" (or some other string) to get whole virus
	
	"""
	
	#make viral chunk
	def __init__(self, viruses, min_len, part = "rand"):

		#get virus to integrate
		self.virus = np.random.choice(list(viruses.keys()))
		#set minimum size for a chunk of virus
		#don't want the size too small (ie only a few base pairs)  
		min_chunk = min_len
		
		if min_chunk>len(viruses[self.virus].seq):
			raise OSError("Viral genome is shorter than the minimum intgration size")
		
		#if we want a random chunk of virus
		if part == "rand":		
			while True:
				self.start = np.random.randint(0, len(viruses[self.virus].seq)-1)
				self.stop = np.random.randint(self.start+1, len(viruses[self.virus].seq))
				if self.stop-self.start>min_chunk:
					break					
		#if we want the whole virus
		else:
			self.start = 0
			self.stop = len(viruses[self.virus].seq)
		if np.random.uniform() > 0.5:
			self.ori = "f" #define orientation
			self.bases = viruses[self.virus][self.start:self.stop] #get bases to insert
		else:
			self.ori = "r" #define orientation
			self.bases = viruses[self.virus][self.start:self.stop].reverse_complement() #get bases to insert remove this one
		#construct dictionary with keys 'bases', 'ori' and 'coords'
		#use to keep track of order if 
		self.pieces = {0:{"bases":self.bases, "ori":self.ori, "coords":(self.start, self.stop)}}
		#store information about rearrangements
		self.isSplit = False #has chunk been split?
		self.isRearranged = False #has chunk been rearranged?
		self.deletion = False #does this chunk have a deletion
		
	def split(self, n):
		#split a part of a virus into n random parts
		#for later rearrangement or deletion
		
		if self.isSplit is True:
			print("Already split")
			return
		self.isSplit = True
		
		#need to check that we have enough bases to leave at least one base in each section
		if self.stop - self.start < n:
			print("Not enough bases to split")
			return
		
		#make orientations of each part the same as the original ori
		oris = [self.pieces[0]["ori"] for i in range(n)]
		
		#get n random coordinates to be breakpoints (plus original start and stop)
		breaks = [self.start, self.stop] + list(np.random.randint(self.start+1, self.stop-1, n-1))
				
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
		
		#don't rearrange chunk if it's already been rearranged
		if self.isRearranged is True:
			raise ValueError("Chunk has already been rearranged!")
			
		self.isRearranged = True #rearranged is now true
		
		#get new order, and make sure it's not the same as before
		order = np.arange(n)
		while True:
			np.random.shuffle(order)
			if not(all([i==j for i, j in zip(order, range(len(order)))])):
				break
		
		#replace keys with shuffled numbers
		self.pieces = {order[key]:value for key, value in self.pieces.items()}
			
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
			
		#don't rearrange chunk if it's already been rearranged
		if self.isRearranged is True:
			raise ValueError("Chunk has already been rearranged!")
			
		self.deletion = True #deletion is now true	
		
	
		#get piece to delete - use middle piece
		key_del = np.around((max(self.pieces.keys()) - min(self.pieces.keys()))/2)
		

		#replace keys with shuffled numbers
		del self.pieces[key_del]
		
		#return the number of deleted fragments
		return 

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

		for i in range(1,len(intCoords)):

			#consider if an integration site was moved due to an overlap 
			if int_list[i].newpoint > -1: 
				i_site = int_list[i].newpoint
			else: 
				i_site = int_list[i].hPos 
				i_site = intCoords[i][0] 
			for j in range(0,i):
				#consider if an integration site was moved due to an overlap 
				if int_list[j].newpoint > -1: 
					j_site = int_list[j].newpoint
				else: 					
						j_site = int_list[j].hPos
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
						new_coord1 = intCoords[j][0]+int_list[i].numBases + shift
						#stop coordinate of the integraton 
						new_coord2 = intCoords[j][1]+int_list[i].numBases + shift  
						intCoords[j] = (new_coord1,new_coord2)
						print("overlap adjust") 
						
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
