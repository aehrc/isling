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
	np.random.seed(1000)

	#types of insertions
	insertion_types = [insertWholeVirus, insertViralPortion, insertWholeRearrange, insertWithDeletion]
	#junction_types = ['gap', 'overlap', 'clean']
	
	host_ints = []
	
	#open output file for writing
	handle = open(args.locs, "w+")
	header = "\t".join(["hChr",
						"hPos",
						"num_fragments",
						"virus",
						"oris",
						"rearranged",
						"vStart",
						"vStart",
						"breakpoints",
						"vBases\n"
							])
	handle.write(header)
	
	insertion_types[0](host, virus, host_ints, handle, sep=args.sep)
	insertion_types[1](host, virus, host_ints, handle, sep=args.sep)
	insertion_types[2](host, virus, host_ints, handle, sep=args.sep)
	
	
	#insertion_types[3](host, virus, host_ints, handle)
	
	handle.close()

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
			return None
	
	#write to output file
	currentInt.writeIntegration(filehandle)
	
	print(currentInt)
	
	return currentInt

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
			return None
	
	#write to output file
	currentInt.writeIntegration(filehandle)
	
	print(currentInt)
	
	return currentInt

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
			return None
	
	#write to output file
	currentInt.writeIntegration(filehandle)
	
	print(currentInt)
	
	return currentInt

def insertWithDeletion(host, viruses, int_list, filehandle, sep=5):
	""" Inserts n portions of viral DNA into host genome"""
	
	#make integration object to store properties of integration
	currentInt = Integration(host)
	
	
	#write to output file
	currentInt.writeIntegration(filehandle)
	
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
		"""
		initalize integration at a random place in the host genome
		only one integration, but may be of a rearranged chunk
		"""
		self.host = host
		self.chr = np.random.choice(list(host.keys()))
		self.hPos = np.random.randint(1, len(host[self.chr].seq))
		self.fragments = 0
		
	def addFragment(self, viruses, part = "rand"):
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
		self.chunk = ViralChunk(viruses, part)
		
		#overlaps and gaps at each end of the integration
		self.overlaps = (0,0) #TODO - allow for gaps and overlaps
		
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
		
		#overlaps and gaps at each end of the integration
		self.overlaps = (0,0) #TODO - allow for gaps and overlaps
			
		#add a rearranged fragment
		self.chunk = ViralChunk(viruses, part)
		self.chunk.rearrange(n)
		
	def doIntegration(self, filehandle):
		"""
		Inserts viral DNA at position index in host DNA
		"""
	
		#do overlap, gap, clean junction
		####TODO####
		pass
		#if self.overlap 
		#return host[:index] + "".join(self.bases) + host[index:]
	
	def writeIntegration(self, filehandle):
		#write information about the current integration to an output file
		#pass in an open filehandle for writing
		
		if self.fragments == 0:
			print("No fragments yet to write")
			return
			
		#get list of oris, breakpoints, bases
		oris = [self.chunk.pieces[i]["ori"] for i in range(self.fragments)]
		coords = [self.chunk.pieces[i]["coords"] for i in range(self.fragments)]
		bases = [self.chunk.pieces[i]["bases"] for i in range(self.fragments)]

		#construct lines to write for current integration
		line = "\t".join([self.chr, 
							str(self.hPos), 
							str(self.fragments), 
							self.chunk.virus,
							"/".join(oris),
							str(self.chunk.isRearranged),
							str(self.chunk.start),
							str(self.chunk.stop),
							"/".join([f"{start}-{stop}" for (start, stop) in coords]),
							"/".join([str(x.seq) for x in bases])+"\n"])
				
		filehandle.write(line)

		#match format to output of pipeline

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
	
		#get types depending on overlaps
		self.types = ["none" if bp == 0 else "gap" if bp < 0 else "overlap" for bp in overlapBPs]
		

	def __str__(self):
		if self.fragments == 0:
			return f"Integration on chromosome {self.chr}, position {self.hPos}"
		elif self.fragments == 1:
			oris = [self.chunk.pieces[i]["ori"] for i in range(self.fragments)]
			coords = [self.chunk.pieces[i]["coords"] for i in range(self.fragments)]
			bases = [self.chunk.pieces[i]["bases"].seq for i in range(self.fragments)]
			return f"Integration of {self.fragments} viral fragument on chromosome {self.chr}, position {self.hPos}.\n" \
				f"Integrated virus: {self.chunk.virus}, virus bases: {bases}, virus orientations: {oris}, coords {coords}"	
		else:
			oris = [self.chunk.pieces[i]["ori"] for i in range(self.fragments)]
			coords = [self.chunk.pieces[i]["coords"] for i in range(self.fragments)]
			bases = [self.chunk.pieces[i]["bases"].seq for i in range(self.fragments)]
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
		
		#if we want a random chunk of virus
		if part == "rand":		
			self.start = np.random.randint(0, len(viruses[self.virus].seq))
			self.stop = np.random.randint(self.start+1, len(viruses[self.virus].seq))
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
		#delete a random piece from a viral chunk of n pieces
		
		pass
	
		print(self.bases)
if __name__ == "__main__":
	main(sys.argv[1:])