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
	#overlaps/gaps at each junction
# reports all integration sites relative to the host genome, independent of other integrations
# so if there are two integrations on the same chromosome, both positions in host genome will be relative to reference


###import libraries
from Bio import SeqIO
from Bio.Alphabet.IUPAC import unambiguous_dna
import argparse
import sys
import os
import random
import numpy as np

### main
def main(argv):

	#get arguments
	parser = argparse.ArgumentParser(description='simulate viral insertions')
	parser.add_argument('--host', help='host fasta file', required = True)
	parser.add_argument('--virus', help = 'virus fasta file', required = True)
	parser.add_argument('--integrated', help = 'output fasta file', required = True)
	parser.add_argument('--locations', help = 'output csv with integration locations', required = True)
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
	random.seed(1000)

	#types of insertions
	insertion_types = [insertWholeVirus, insertViralPortion, insertRearrange, insertRearrangeDeletion]
	junction_types = ['gap', 'overlap', 'clean']
	
	host_int = []
	
	insertion_types[0](host, host_int, virus)
	insertion_types[1](host, host_int, virus)
	insertion_types[2](host, host_int, virus)
	insertion_types[3](host, host_int, virus)

def insertWholeVirus(host, host_int, viruses):
	"""Inserts whole viral genome into host genome"""
	#make integration object to store properties of integration
	currentInt = Integration(host)
	
	#get one viral chunk
	currentInt.addFragment(viruses, part = "whole")
	
	#do integration
	print(currentInt)
	
	return currentInt, host_int

def insertViralPortion(host, host_int, viruses):
	"""Inserts portion of viral DNA into host genome"""
	
	#make integration object to store properties of integration
	currentInt = Integration(host)
	
	#get one viral chunk
	currentInt.addFragment(viruses, part = "rand")
	
	#do integration

	
	print(currentInt)

	return currentInt, host_int

def insertRearrange(host, host_int, viruses, n = 2):
	""" Inserts a single portion of viral DNA with n rearrangements """
	
	#make integration object to store properties of integration
	currentInt = Integration(host)
	
	#get one viral chunk
	currentInt.addFragment(viruses)
	
	#divide up chunk into n fragments
	
	return currentInt, host_int

def insertRearrangeDeletion(host, host_int, viruses, n = 2):
	""" Inserts n portions of viral DNA into host genome"""
	
	#make integration object to store properties of integration
	currentInt = Integration(host)
	
	for i in range(1, n):
		currentInt.addFragment(viruses)
	
	return currentInt, host_int

def doIntegration(host, index, virus):
	"""Inserts viral DNA at position index in host DNA"""
	
	#do overlap, gap, clean junction
	
	return host[:index] + virus + host[index:]
	
	
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
	
class Integration:
	""" 
	Class to store the properties of integrations 
	"""
	
	def __init__(self, host):
		#initialise with random place in the host genome
		self.host = host
		self.chr = random.choice(list(host.keys()))
		self.hStart = random.randint(1, len(host[self.chr].seq))
		self.hStop = self.hStart #to store location of integration stop
		
		#initialise with no fragments
		self.fragments = 0 #to store the number of fragments of virus integrated
		self.viruses = [] #to store the viruses that have been integrated
		self.vParts = [] #to store the start and stop of each integrated viral fragment
		self.oris = [] #to store orientations of each fragment
		self.bases = [] #to store bases inserted
		
	def addFragment(self, viruses, part = "rand"):
		#add a simple fragment
		#get viral chunk
		chunk = ViralChunk(viruses, part)
		
		#check that chunk doesn't overlap with any current chunks
		
		#add viral chunk
		self.fragments += 1 #increment number of fragments
		self.viruses.append(chunk.virus)
		self.vParts.append([chunk.start, chunk.stop])
		self.oris.append(chunk.ori)
		self.bases.append(chunk.bases)
		
	def addRearrange(self, viruses, part = "rand"):
		
		#add a rearranged fragment
		chunk = ViralChunk(viruses, part)
		
		#get number of fragments to rearrange into
		#draw from a gamma distribution with k = 2 and theta = 2
		
		n = np.random.gamma(2,2)
		
		chunk.rearrange(n)
		
		
	def __str__(self):
		if self.fragments == 0:
			return "Integration on chromosome {}, position {}".format(self.chr, self.hStart)
		else:
			return "Integration of {} viral fraguments on chromosome {}, position {}.  Integrated viruses: {}, virus bases: {}, virus orientations: {}".format(self.fragments, self.chr, self.hStart, self.viruses, self.vParts, self.oris)
		
class ViralChunk:
	"""
	Class with a chunk of virus
	Default is to get a random chunk (part = 'rand')
	use part = "whole" (or some other string) to get whole virus
	
	"""
	
	#make viral chunk
	def __init__(self, viruses, part = "rand"):
	
		#get virus to integrate
		self.virus = random.choice(list(viruses.keys()))
		
		#if we want a random chunk of virus
		if part == "rand":		
			self.start = random.randint(0, len(viruses[self.virus].seq))
			self.stop = random.randint(self.start, len(viruses[self.virus].seq))
		#if we want the whole virus
		else:
			self.start = 1
			self.stop = len(viruses[self.virus].seq)
	
		#forward or reverse?
		if random.random() > 0.5:
			self.ori = "f" #define orientation
			self.bases = viruses[self.virus].seq[self.start:self.stop] #get bases to insert
		else:
			self.ori = "r" #define orientation
			self.bases = viruses[self.virus].seq[self.start:self.stop].reverse_complement() #get bases to insert

		
		#store information about rearrangements
		self.rearranged = False #has chunk been rearranged?
		self.breakpoints = [] #store breakpoints
		self.oris = [self.ori] #store orientations of rearranged parts of chunk
			
	def rearrange(self, n):
		#rearrange a chunk in n parts
		#pick n random breakpoints, and extract sequence of each
		#recombine in a random order, with equal probability of each
		#part being in forward or reverse order
		
		#don't rearrange chunk if it's already been rearranged
		if self.rearranged is True:
			raise ValueError("Chunk has already been rearranged!")
		
		self.rearranged = True #rearranged is now true
		
		#get n random coordinates to be breakpoints
		#and n random orientations
		while len(self.breakpoints) < (n - 1):
			self.breakpoints.append(random.randint(self.start, self.stop))
			#get random number for orientations
			if random.random() > 0.5:
				self.oris.append("f")
			else:
				self.oris.append("r")
		
		#add start and stop to list of breakpoints to make getting sequence easier
		self.breakpoints.append(self.start)
		self.breakpoints.append(self.stop)
				
		#sort breakpoints to make them easier to work with
		self.breakpoints.sort()
		
		#rearrange chunk sequence
		temp = self.bases[0:0]
		for i in range(n):
			#get start and stop for this breakpoint
			b_start = self.breakpoints[i]
			b_stop = self.breakpoints[i+1]
			
			#add sequence to temp
			if self.oris[i] == 'f':
				temp += viruses[self.virus].seq[b_start:b_stop]
			else:
				temp += viruses[self.virus].seq[b_start:b_stop].reverse_complement()
				
		self.bases = temp
	

if __name__ == "__main__":
	main(sys.argv[1:])