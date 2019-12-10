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
import pandas as pd
import argparse
import sys
import os
import numpy as np

###
max_attempts = 5 #maximum number of times to try to place an integration site

### main
def main(argv):

	#get arguments
	parser = argparse.ArgumentParser(description='simulate viral insertions')
	parser.add_argument('--host', help='host fasta file', required = True)
	parser.add_argument('--virus', help = 'virus fasta file', required = True)
	parser.add_argument('--ints', help = 'output fasta file', required = True)
	parser.add_argument('--locs', help = 'output csv with integration locations', required = True)
	parser.add_argument('--fasta', help = 'output fasta of integrated host genome', required = False)
	parser.add_argument('--sep', help = 'integrations must be seperated by this many bases', required=False, default=5)
	parser.add_argument('--min_len', help = 'minimum length of integerations', required=False, default=5)
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
	np.random.seed(500)

	#types of insertions
	insertion_types = [insertWholeVirus, insertViralPortion, insertWholeRearrange, insertWithDeletion]
	
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
	int_num = 2
	
	#intialise how many episomal sequences included in the outputted fasta file 
	epi_num = 5
	
	#intialise after how many intergrations the number of integrations performed is reported to the user 
	int_report = 100
	
	
	print("\nNUMBER OF INTEGRATIONS TO INSERT: "+str(int_num))
	
	#integration loop 
	for i in range(0,int_num):
		rand_int =  np.random.randint(0,4)
		host_ints, host_fasta = insertion_types[rand_int](host_fasta, virus, host_ints, handle, sep=args.sep)
		if i % int_report == 0 and i != 0: 
			print(str(i) +" integrations complete...")
			
	print("\nNUMBER OF EPISOMES: "+str(epi_num))
	
	#adds viral sequences as 'episomes'  
	virus_key = list(virus.keys())

	for i in range(0,epi_num): 
		#host_fasta.append(virus.get('virus'))
		host_fasta['episome '+str(i)] = virus.get(virus_key[0])
			
	print("\n***INTEGRATIONS COMPLETE***")
	print(host_fasta)
	handle.close()
	
	#save statistics on the integration 
	stats = Statistics.saveStats(host_ints)
	stats.to_csv("host_insertions.csv",sep = '\t') 
	
	#save integrated host sequence 
	with open(args.ints, 'w') as handle: 
    		SeqIO.write(host_fasta.values(), handle, 'fasta')
    		handle.close()
    		print("\nIntegrated host saved as "+str(args.ints))
    		print("Details of sequences integrated saved as "+str(args.locs))
    		print("Details of where intgrations lie in host sequence saved as "+"host_insertions.csv")


def insertWholeVirus(host, viruses, int_list, filehandle, sep=5):
	"""Inserts whole viral genome into host genome"""
	
	#get positions of all current integrations
	currentPos = [int.hPos for int in int_list]
	#make sure that new integration is not within args.sep bases of any current integrations
	attempts = 0
	while True:
		#get one viral chunk
		currentInt = Integration(host)
		currentInt.addFragment(viruses, part = "whole")
		attempts += 1
		
		#check that all integrations are sep away from new integration
		if all([ abs(currentInt.hPos-i) > sep for i in currentPos ]):
			break
		#if we've made a lot of attempts, give up
		elif attempts > max_attempts:
			return int_list, host
	
	#do integration
	host = currentInt.doIntegration(host, int_list)
	
	#write to output file
	currentInt.writeIntegration(filehandle)
	
	#append to int_list
	int_list.append(currentInt)
	
	return int_list, host

def insertViralPortion(host, viruses, int_list, filehandle, sep=5):
	"""Inserts portion of viral DNA into host genome"""
	
	#get positions of all current integrations
	currentPos = [int.hPos for int in int_list]
	
	#make sure that new integration is not within args.sep bases of any current integrations
	attempts = 0
	while True:
		#get one viral chunk
		
		currentInt = Integration(host)
		currentInt.addFragment(viruses, part = "rand")
		attempts += 1
		
		#check that all integrations are sep away from new integration
		if all([ abs(currentInt.hPos-i) > sep for i in currentPos ]):
			break
		#if we've made a lot of attempts, give up
		elif attempts > max_attempts:
			return int_list, host
			
	#do integration
	host = currentInt.doIntegration(host, int_list)
	
	#write to output file
	currentInt.writeIntegration(filehandle)
	
	#append to int_list
	int_list.append(currentInt)
	
	return int_list, host

def insertWholeRearrange(host, viruses, int_list, filehandle, sep=5):
	""" Inserts a single portion of viral DNA with n rearrangements """
	
	#get positions of all current integrations
	currentPos = [int.hPos for int in int_list]
	
	#make sure that new integration is not within args.sep bases of any current integrations
	attempts = 0
	while True:
		#get one viral chunk
		currentInt = Integration(host)
		currentInt.addRearrange(viruses, part = "whole")
		attempts += 1
		
		#check that all integrations are args.sep away from new integration
		if all([ abs(currentInt.hPos-i) > sep for i in currentPos ]):
			break
		#if we've made a lot of attempts, give up
		elif attempts > max_attempts:
			return int_list, host
			
	#do integration
	host = currentInt.doIntegration(host, int_list)
	
	#write to output file
	currentInt.writeIntegration(filehandle)
	
	#append to int_list
	int_list.append(currentInt)
	
	return int_list, host

def insertWithDeletion(host, viruses, int_list, filehandle, sep=5):
	""" Inserts n portions of viral DNA into host genome"""
	
	#get positions of all current integrations
	currentPos = [int.hPos for int in int_list]
	
	#make sure that new integration is not within args.sep bases of any current integrations
	attempts = 0
	while True:
		#get one viral chunk
		currentInt = Integration(host)
		currentInt.addDeletion(viruses, part = "whole")
		attempts += 1
		
		#check that all integrations are args.sep away from new integration
		if all([ abs(currentInt.hPos-i) > sep for i in currentPos ]):
			break
		#if we've made a lot of attempts, give up
		elif attempts > max_attempts:
			return int_list, host
	
	#do integration
	host = currentInt.doIntegration(host,int_list)
	
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
 		 
	def addFragment(self, viruses, part = "rand"):
		"""
		add a viral fragment to this integration
		"""
		
		#check there isn't alreafbdy a fragmentfg
		if self.fragments > 0:
			print("Fragment has already been added!")
			return
			
		#add a simple fragment
		self.fragments = 1
		
		#get viral chunk
		self.chunk = ViralChunk(viruses, part)

	def addRearrange(self, viruses, part = "rand"):
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
		self.chunk = ViralChunk(viruses, part)
		self.chunk.rearrange(n)
	
	def addDeletion(self, viruses, part = "rand"):
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
		self.chunk = ViralChunk(viruses, part)
		self.chunk.delete(n)
		self.fragments -= 1


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
		
	def overlapPoint(self,left_site,right_site):
		"""
		Finds point where there is homology between the human and viral sequences 
		"""
		overlap_point = 0
		#if no homologous region is found we skip the integration 
		if left_site ==-1 and right_site == -1: 
			self.chunk.bases = ""
			self.inserted = False 
		
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
		return overlap_point  

	def createLeftOverlap(self,host,int_list):
		"""
		Handles overlaps on the left of a viral chunk. Left and right are treated as different functions as different operations must be performed. 
		Works by finding closest regions of homology on left side of the randomly selected integration point. It checks if the homologous region is caused by an existing integration and concatenates the search range and attempts to find a homologous region again if the homology was caused by an existing integration. This is repeated for the right side of te integration point. The region closest to randomly selected integration point is then used to insert the viral chunk. 
		"""
 		
		viral_chunk = self.chunk.bases

		l_overlap = str(viral_chunk[:-self.overlaps[0]].seq) #sequence left of the integration site 
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
				left_seq = host[self.chr][:left_search]
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
		
		overlap_point = self.overlapPoint(left_site,right_site)
		int_start = overlap_point
		int_stop = overlap_point 

		if left_site != -1 and right_site != -1: 
		#only if an integration occurs does the overlapping region get removed 
			int_start =  overlap_point + self.overlaps[0] #TODO recheck this 
		return int_start,int_stop 
		
	def createRightOverlap(self,host,int_list): 
		"""
		Handles overlaps on the right of a viral chunk. Same as above with operations applicable to right end 
		""" 
		viral_chunk = self.chunk.bases
		dont_integrate = self.dontIntegrate(int_list)
		
		r_overlap = str(viral_chunk[len(viral_chunk)+self.overlaps[1]:].seq) #sequence right of the integration site 
		left_site = -1
		right_site = -1
	
		#find homologous region on the left side 
		left_seq = host[self.chr][:self.hPos].seq

		while left_site<0:
			left_search = str(left_seq).rfind(r_overlap)
			if left_search == -1:
				break 
			if left_search in dont_integrate: # check homology is not from a previous integration
				left_seq = host[self.chr][:left_search]	
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
				right_seq = host[self.chr][right_search:]
			else: 
				right_site = right_search
		
		overlap_point = self.overlapPoint(left_site,right_site)
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
		
		
			
	def doIntegration(self, host, int_list):
		"""
		Inserts viral DNA (self.bases) at position self.hPos in host[self.chr]
		"""
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
			(int_start,int_stop) = self.createLeftOverlap(host,int_list)
			 
		#adjust sequences with right overlap  		
		if self.overlaps[1]<0:
			(int_start,int_stop) = self.createRightOverlap(host,int_list)
		
		host[self.chr] = host[self.chr][:int_start] + \
		 				"".join(self.chunk.bases) + \
		 				host[self.chr][int_stop:]
		
		#set the starting and stopping points of the performed integration 
		self.setStopStart(int_start,int_stop) 
		 				
		return host


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
		""" provides coordinates of where viral DNA starts and ends"""

		
		if self.overlaps[0]>0: #ie there is a gap 
			self.start = int_start+self.overlaps[0]
		elif self.overlaps[0]<0: #ie there is an overlap
			self.start = int_start-self.overlaps[0]
		else: 
			self.start = int_start 
			
		if self.overlaps[1]>0: #ie there is a gap
			self.stop = int_stop-self.overlaps[1]
		elif self.overlaps[1]<0: #ie there is an overlap 
			self.stop = int_stop+self.overlaps[1]
		else: 
			self.stop = int_stop 
			
		self.stop = self.stop+len(self.chunk.bases) 

		
	def getStartStop(self):
		"""
		return starts and stops of each end of viral integration, 
		accounting for any gaps or overlaps
		two integration site for each integrated bit of virus
		"""
		
		#human coordinates - 0-based
		self.hStarts = [self.hPos, self.hPos+self.overlaps[1]]
		self.hStops = [self.hPos+1+self.overlaps[0], self.hPos+1]
		
		self.vStarts = []
		self.vStops = []
		#viral coordinates depend on orientation
		if self.oris[0] == "f":
			self.vStarts.append(self.vParts[0][0])
		
		
		#vStart and stop are relative to virus
		#so if forward, vStart is first site relative to host
		#otherwise is second
	
		if ori == "f":
			coords1 = [pos, pos+overlapBPs[0]+1].sort()
			coords2 = [pos, pos+overlapBPs[1]+1].sort()
		else:
			coords1 = [pos, pos+overlapBPs[0]+1].sort()
			coords2 = [pos, pos+overlapBPs[1]+1].sort()
	
		#rearrange starts and stops for output
		self.starts = (coords1[0], coords2[0])
		self.stops = (coords1[1], coords2[1])
	
		#orientations are always human-virus and virus-human
		self.oris = ["hv", "vh"]
				
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
	def __init__(self, viruses, part = "rand"):
	
		#get virus to integrate
		self.virus = np.random.choice(list(viruses.keys()))

		#set minimum size for a chunk of virus
		#don't want the size too small (ie only a few base pairs)  
		min_chunk = 50 
		
		#if we want a random chunk of virus
		if part == "rand":		
			while True: 
				self.start = np.random.randint(0, len(viruses[self.virus].seq))
				if self.start+1 != len(viruses[self.virus].seq):
					break
					
			#ensure size of viral chunk is greater than the minimum size 
			while True: 
				self.stop = np.random.randint(self.start+1, len(viruses[self.virus].seq))
				if self.stop-self.start>min_chunk:
					break
					
				
		#if we want the whole virus
		else:
			self.start = 0
			self.stop = len(viruses[self.virus].seq)
	
		#forward or reverse?
		if np.random.uniform() > 0.5:
			self.ori = "f" #define orientation
			self.bases = viruses[self.virus][self.start:self.stop] #get bases to insert
		else:
			self.ori = "r" #define orientation
			self.bases = viruses[self.virus][self.start:self.stop].reverse_complement() #get bases to insert

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
		
		intCoords = [(int.start,int.stop) for int in int_list]
		for i in range(1,len(intCoords)):
			for j in range(0,i):
				if int_list[i].hPos < int_list[j].hPos: 
					#intCoords[j] = intCoords[j]+int_list[i].numBases
					 new_coord1 = intCoords[j][0]+int_list[i].numBases
					 new_coord2 = intCoords[j][1]+int_list[i].numBases
					 intCoords[j] = (new_coord1,new_coord2) 
					
		return intCoords
		
	def saveStats(int_list):
		"""Save a csv file with the coordinates of the adjusted insertions so when we create artifical reads we can predict which reads contain viral DNA""" 
		num_ints = len(int_list) 
		
		int_sites = Statistics.intList(int_list[1-num_ints],int_list)
		int_coords = Statistics.adjustedStopStart(int_list[1-num_ints],int_list) 
		
		int_start = []
		int_stop = []
		
		for i in range(0,num_ints): 
			c1,c2 = int_coords[i]
			int_start.append(c1)
			int_stop.append(c2)
		
		stats_df = pd.DataFrame({"Integration sites":int_sites,"Start point":int_start,"Stop point": int_stop}) 
		
		return stats_df 
		 

if __name__ == "__main__":
	main(sys.argv[1:])
