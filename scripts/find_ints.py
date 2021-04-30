#!/usr/bin/env python3

import argparse
import pysam
import csv
import re
import array
import itertools
import pdb


# note python 3.7+ is required

default_header = ['Chr', 'IntStart', 'IntStop', 'VirusRef', 'VirusStart', 'VirusStop', 
		'NoAmbiguousBases', 'OverlapType', 'Orientation', 'VirusOrientation', 
		'HostSeq', 'VirusSeq', 'AmbiguousSeq', 'HostEditDist', 'ViralEditDist', 
		'TotalEditDist', 'PossibleHostTranslocation', 'PossibleVectorRearrangement', 
		'HostAmbiguousLocation', 'ViralAmbiguousLocation', 'Type', 
		'HostMapQ', 'ViralMapQ', 'ReadID', 'ReadSeq', 'AlternativeInts']

def main():
	
	parser = argparse.ArgumentParser(description='Identify integrations by comparing a host and viral alignment')
	parser.add_argument('--host', help="host alignment bam (query-sorted)", required=True)
	parser.add_argument('--virus', help="virus alignment bam (query-sorted)", required=True)
	parser.add_argument('--integrations', help='file to output integrations to', default='integrations.tsv')
	parser.add_argument('--map-thresh', help="threshold for length of mapped region", default=20, type=int)
	parser.add_argument('--mean-template-length', help="mean template length for library", default=0, type=int)
	parser.add_argument('--tolerance', help="tolerance for short CIGAR operations (MIS ops with a combined length equal to or shorter than this will be combined into the nearest mapped region)", type=int, default=3)
	parser.add_argument('--verbose', '-v', help="print extra messages?", action='store_true')
	args = parser.parse_args()

	pair = AlignmentFilePair(args.host, args.virus, args.map_thresh, 
								args.mean_template_length, args.tolerance, args.verbose)
	
	pair.find_integrations()
	
	pair.write_integrations(args.integrations)
		
	pair.close()
	
		
class AlignmentFilePair:
	""" A class to hold host and viral alignments, and look for integrations """
	
	def __init__(self, host, virus, map_thresh, tlen, tol = 3, verbose = False):
		""" Open host and viral alignments, and check for query sorting """
		
		self.host = AlignmentFile(host, verbose)
		self.virus = AlignmentFile(virus, verbose)	
		
		self.href_lens = self.host.get_reference_lengths()
		self.vref_lens = self.virus.get_reference_lengths()
		self.map_thresh = map_thresh
		self.verbose = verbose
		self.tlen = tlen
		self.tol = tol
		self.ints = []
		
		self.default_header = default_header
		
		if self.tlen == 0:
			print(f"warning: template length is zero - the position for discordant integraiton sites will be at the end of the mapped read")
		
		if self.tlen < 0:
			raise ValueError("Template length must be a positive integer")
		
		if self.tol < 0:
			raise ValueError("Tolerance must be a positive integer")
		
		if self.verbose:
			print(f"using host bam {host} and virus bam {virus}")
			

	def find_integrations(self):
		""" Check for integrations """
		
		counter = 0
		progress = 100000
		
		# iterate over reads that have alignments to both host and virus
		for alns in self._get_aligns():
		
			counter += 1
			if counter % progress == 0:
				print(f"checked {counter} reads")
			
			if self.curr_host is None or self.curr_virus is None:
				break
				
			if self.verbose:
				print(f"checking {self.curr_host[0].qname}")
				
			# collect read 1 alignments
			host_r1_primary  = self.curr_host.get_primary_r1()
			host_r1_not_primary = self.curr_host.get_not_primary().get_read1()
			
			virus_r1_primary =  self.curr_virus.get_primary_r1()
			virus_r1_not_primary = self.curr_virus.get_not_primary().get_read1()

			# check if primary read 1 alignments appear to be chimeric	
			if host_r1_primary is not None and virus_r1_primary is not None:
				integration = self._is_chimeric(host_r1_primary, virus_r1_primary,
												host_r1_not_primary, virus_r1_not_primary)
												
				if integration is not None:
					self.ints.append(integration)
					integration.get_properties()

			# collect read2 alignments
			host_r2_primary = self.curr_host.get_primary_r2()
			host_r2_not_primary = self.curr_host.get_not_primary().get_read2()
			
			virus_r2_primary = self.curr_virus.get_primary_r2()
			virus_r2_not_primary = self.curr_virus.get_not_primary().get_read2()			
			
			# check if primary read 1 alignments appear to be chimeric	
			if host_r2_primary is not None and virus_r2_primary is not None:
				integration = self._is_chimeric(host_r2_primary, virus_r2_primary,
												host_r2_not_primary, virus_r2_not_primary)
				if integration is not None:
					self.ints.append(integration)
					integration.get_properties()						
						
			# if there are alignments that are neither read 1 nor read 2, check if chimeric
			host_not12_primary = self.curr_host.get_primary_not_r1_or_r2()
			host_not12_not_primary = self.curr_host.get_not_primary().get_not_read1_or_read2()
			
			virus_not12_primary = self.curr_virus.get_primary_not_r1_or_r2()
			virus_not12_not_primary = self.curr_virus.get_not_primary().get_not_read1_or_read2()
			
			if host_not12_primary is not None and virus_not12_primary is not None:
				integration = self._is_chimeric(host_not12_primary, virus_not12_primary,
												host_not12_not_primary, virus_not12_not_primary)
				if integration is not None:
					self.ints.append(integration)
					integration.get_properties()
					
			# check for a discordant pair
			if all([i is not None for i in [host_r1_primary, host_r2_primary, 
											virus_r1_primary, virus_r2_primary]]):
				# check if integration
				integration = self._is_discordant(host_r1_primary, host_r2_primary, 
											virus_r1_primary, virus_r2_primary,
											host_r1_not_primary, host_r2_not_primary,
											virus_r1_not_primary, virus_r2_not_primary)
				if integration is not None:
					self.ints.append(integration)
					integration.get_properties()
		
		print(f"found {len(self.ints)} integrations")
		
	def write_integrations(self, filename, header = None):
		""" Write integrations to file """
		
		if header is None:
			header = self.default_header
			
		try:
			with open(filename, 'w') as filehandle:
				out = csv.DictWriter(filehandle, fieldnames = header, delimiter='\t')
				out.writeheader()
				for integration in self.ints:
					out.writerow(integration.get_properties())
		except AttributeError:
			print("No file written: integrations have not been detected yet!")
			return
			
		print(f"saved output to {filename}")
				
				
	def _is_chimeric(self, hread, vread, hsec, vsec):	
		""" 
		Takes two alignments for the same read, one from the host and one from the virus, 
		and decides if the read appears to be a simple chimera.  A 'simple chimera' means
		that one part of the read is accounted for by the host alignment, and another
		pairt is accounted for by the viral alignment,
		for example (130M20S for host and 20M130S for virus).  There might be some overlap
		or gap between the host and viral alignments, for example 132M18S for host and
		22M132S for virus would indicate a 2 bp overlap
		"""
		try:
			return ChimericIntegration(hread, vread, 
										self.host.aln.header, 
										self.virus.aln.header,
										hsec, vsec,
										self.tol, self.map_thresh)
		except AssertionError:
			return None
			
	def _is_discordant(self, hread1, hread2, vread1, vread2, 
						host_sec1, host_sec2, virus_sec1, virus_sec2,):
		"""
		Checks if a pair looks like an discordant integration
		"""
		
		try:
			return DiscordantIntegration(hread1, hread2, vread1, vread2, 
										self.tol,
										self.host.aln.header, 
										self.virus.aln.header,
										host_sec1, host_sec2, virus_sec1, virus_sec2,
										self.map_thresh, self.tlen
										)
		except AssertionError:
			return None

	def _get_aligns(self):
		""" A generator to get alignments for the same read from both host and virus"""

		# collect alignments from file
		self.curr_host = self.host.collect_next_read_alignments()
		self.curr_virus = self.virus.collect_next_read_alignments()

		# continue until we've either found the same read or reached the end of one of the files
		while True:
		
			if self.curr_host is None or self.curr_virus is None:
				return
		
			# check if the query names are the same for host and viral alignments
			if self.curr_host[0].qname == self.curr_virus[0].qname:
				yield
				self.curr_host = self.host.collect_next_read_alignments()
				self.curr_virus = self.virus.collect_next_read_alignments()

			# host qname less than viral qname
			elif self.curr_host[0].qname > self.curr_virus[0].qname:
				if self.verbose:
					print(f"warning: host bam has no alignment for read {self.curr_virus[0].qname} in virus bam, skipping")
				self.curr_virus = self.virus.collect_next_read_alignments()

			# viral qname less than host qname
			elif self.curr_virus[0].qname > self.curr_host[0].qname:
				if self.verbose:
					print(f"warning: virus bam has no alignments for read {self.curr_host[0].qname} in host bam, skipping")
				self.curr_host = self.host.collect_next_read_alignments()

	def close(self):
		""" Close both alignments """

		self.host.close()
		self.virus.close()	

class AlignmentFile:
	""" Basically just pysam.AlignmentFile, but with a few extra methods """
	
	def __init__(self, bam_path, verbose):
		""" Open file with pysam, check if query-sorted """
		
		# open file
		self.filename = bam_path
		self.aln = pysam.AlignmentFile(bam_path,'rb')
		self.end = False
		self.verbose = verbose

		# check for query sort 
		error = f"The bam file {bam_path} does not appear to be query-sorted.  Please sort (i.e. samtools sort -n -o <sorted> <unsorted>) and try again"
		if 'HD' not in self.aln.header.to_dict():
			raise ValueError(error)
		if 'SO' not in self.aln.header.to_dict()['HD']:
			raise ValueError(error)			
		if self.aln.header.to_dict()['HD']['SO'] != 'queryname':
			raise ValueError(error)

		# get first read from file
		try:
			self.curr = next(self.aln)
		except StopIteration:
			self.end = True
			self.curr = None


	def collect_next_read_alignments(self):
		""" Collect the alignments (primary, secondary, supplementary) for the next read (pair) in the file """	
		
		alns = AlignmentPool()
			
		qname = self.curr.qname
		
		if self.end:
			return None
		
		# keep getting elements until we either reach a different read or the end of the list
		while (self.curr.qname == qname) and (self.curr.qname is not None):
			# includes R1, R2, primary, secondary, supplementary alignments
			alns.append(self.curr)		
			self._get_next_align()
			if self.end:
				break
				
		return alns

		
	def close(self):
		""" Close alignment file """
		
		self.aln.close()
		
		if self.verbose:
			print(f"closing alignment file {self.filename}")
	
	def get_reference_lengths(self):
		""" Get a dict of name:length for all references """
		
		return {ref:aln for ref, aln in zip(self.aln.references, self.aln.lengths)}
		
	def _get_next_align(self):
		""" Get next alignment from file and assign to self.curr """
		if not self.end:
			try:
				self.curr = next(self.aln)
			except StopIteration:
				self.end = True
				self.curr.qname = None
				if self.verbose:
					print(f"reached end of file {self.filename}")
					
	def _process_align(self):
		"""
		Process alignment to simplify
		
		- change any hard-clipped parts of CIGAR to soft-clip
		"""
		
		for i in range(len(self.curr.cigartuples)):
			if self.curr.cigartuples[i][0] == 5:
				self.curr.cigartuples[i] = (4, self.curr.cigartuples[i][1])

class AlignmentPool(list):
	""" 
	Basically just a list of pysam.AlignedSegment objects, but with a few extra methods.
	
	Should only contain alignments from the same read-pair (i.e. qname should be the same for all reads in pool)
	"""
	
	def add_read_from_XA(self, XA, primary, header):
		""" Add a read based on an XA string (chr,pos,CIGAR,NM) """
		
		# get info from XA string
		XA = XA.split(",")
		assert len(XA) == 4
		
		chr = XA[0]
		ori = XA[1][0]
		pos = int(XA[1][1:]) - 1
		cigar = XA[2]
		nm = int(XA[3])
		
		new_seg = self._create_new_segment_from_primary(primary, header, chr, pos, ori,
														cigar, nm)
		self.append(new_seg)
	
	def add_read_from_SA(self, SA, primary, header):
		""" Add a read based on an SA string (rname,pos,strand,CIGAR,mapQ,NM) """

		# get info from SA string
		SA = SA.split(",")
		assert len(SA) == 6
		
		chr = SA[0]
		pos = int(SA[1]) - 1
		ori = SA[2]
		cigar = SA[3]
		mapq = int(SA[4])
		nm = int(SA[5])
		
		new_seg = self._create_new_segment_from_primary(primary, header, chr, pos, ori,
														cigar, nm, mapq)
		self.append(new_seg)
	
	def get_primary(self):
		""" Return an AlignmentPool with only the primary alignment(s) from an AlignmentPool """
		return AlignmentPool([read for read in self if not read.is_supplementary and not read.is_secondary])
	
	def get_not_primary(self):
		""" Return an AlignmentPool with only the non-primary alignment(s) from an AlignmentPool """
		return AlignmentPool([read for read in self if read.is_supplementary or read.is_secondary])		
		
	def get_read1(self):
		""" Return an AlignmentPool with only read 1 alignment(s) from an AlignmentPool"""
		return AlignmentPool([read for read in self if read.is_read1])
		
	def get_read2(self):
		""" Return an AlignmentPool with only the read 2 alignment(s) from an AlignmentPool """
		return AlignmentPool([read for read in self if read.is_read2])
		
	def get_not_read1_or_read2(self):
		""" Return an AlignmentPool with only the reads that are neither read1 nor read2 from an AlignmentPool """
		return AlignmentPool([read for read in self if not read.is_read1 and not read.is_read2])
		
	def get_primary_r1(self):
		""" Return the primary read 1 alignment as a pysam.AlignedSegment object, or None if there isn't one """
		
		primary_r1 = self.get_primary().get_read1()
		if len(primary_r1) == 0:
			return None
		elif len(primary_r1) == 1:
			return primary_r1[0]	
		else:
			raise ValueError(f"found more than one primary alignment for {self[0].qname} read 1")	
		
	def get_primary_r2(self):
		""" Return the primary read 2 alignment as a pysam.AlignedSegment object, or None if there isn't one """
		
		primary_r2 = self.get_primary().get_read2()
		if len(primary_r2) == 0:
			return None
		elif len(primary_r2) == 1:
			return primary_r2[0]	
		else:
			raise ValueError(f"found more than one primary alignment for {self[0].qname} read 2")	

	def get_primary_not_r1_or_r2(self):
		""" Return the primary read 1 alignment as a pysam.AlignedSegment object, or None if there isn't one """
		
		primary = self.get_primary().get_not_read1_or_read2()
		if len(primary) == 0:
			return None
		elif len(primary) == 1:
			return primary[0]	
		else:
			raise ValueError(f"found more than one primary alignment for {self[0].qname} (not read 1 or read 2)")	
			
	def get_rearrangement_nm(self):
		""" 
		For an AlignmentPool consisting of different alignments of the same read, get an 
		edit distance for the case that those alignments cover the whole read
		
		This is done by collecting non-nested alignments using self.remove_query_nested(),
		and then taking the sum of all the edit distances for those alignments, plus the
		number of bases in the read not accounted for by any alignment (i.e. the number of
		bases in the gaps between those alignments) 
		
		Return a tuple of the edit distance and the non-nested AlignmentPool (since it might
		be different to self after removing nested alignments)
		"""
		
		# check that we're dealing with the same query, and
		# check that all are read1, read2, or neither
		if len(self) > 1:
			assert all([read.query_name == self[0].query_name for read in self[1:]])
			assert all([read.is_read1 == self[0].is_read1 for read in self[1:]])
			assert all([read.is_read2 == self[0].is_read2 for read in self[1:]])
		
		alns = AlignmentPool(self)
		alns.remove_query_nested()
		
		# if there's a gap between the first alignment and the start of the read
		nm = alns[0].query_alignment_start
			
		curr_pos = alns[0].query_alignment_end
		for read in alns:
			
			# if there's a gap between this alignment and the previous one
			if read.query_alignment_start > curr_pos:
				nm += read.query_alignment_start - curr_pos
				
			# get an edit distance for the part of the read
			try:
				nm += read.get_tag('NM')
			except KeyError:
				pass
				
			curr_pos = read.query_alignment_end
				
		# if there's a gap between the last alignment and the end of the read
		nm += alns[-1].query_length - alns[-1].query_alignment_end
		
		return nm, alns
			
	def remove_redundant_alignments(self):
		""" 
		Remove redundant alignments (those that have same pos, CIGAR and strand) 
		Note that two alignments are also redundant if they have the same CIGAR, except that
		soft-clipped operations are replaced by hard-clipped operations
		"""
		
		if len(self) > 1:
			assert all([read.query_name == self[0].query_name for read in self[1:]])
		
		# first, find groups of integrations that are the same
		same = [] # to store pairs of alignments that are redundant
		redundant = set()
		hard_clips = [] # store which alignments contain hard-clips
		for i in range(len(self)-1):
			for j in range(i+1, len(self)):
				if i == j:
					continue
				
				# check for same mapping position
				if self[i].reference_start != self[j].reference_start:
					continue
					
				# check for same strand
				if self[i].is_reverse != self[j].is_reverse:
					continue
				
				# check for same cigar
				i_cigartuples = self._convert_hard_clips_to_soft(self[i].cigartuples)
				j_cigartuples = self._convert_hard_clips_to_soft(self[j].cigartuples)
				if i_cigartuples != j_cigartuples:
					continue
					
				# see if alignments contain hard clips
				if i_cigartuples != self[i].cigartuples:
					hard_clips.append(i)
				if j_cigartuples != self[j].cigartuples:
					hard_clips.append(j)
				
				# check if we've already decided if either of these alignments are
				# redundant with another alignment
				
				if i in redundant:
					pdb.set_trace()
					i_found = [ind for ind in range(len(same)) if i in same[ind]]
					assert len(i_found) == 1
					ind = i_found[0]
					if j not in same[ind]:
						same[ind] = [*same[ind], j]
				elif j in redundant:
					pdb.set_trace()
					j_found = [ind for ind in range(len(same)) if i in same[ind]]
					assert len(j_found) == 1
					if i not in same[ind]:
						same[ind] = [*same[ind], i]				
				else:
					same.append([i, j])
				
				redundant.add(i)
				redundant.add(j)
				
		# remove integrations from groups of redundant alignments, leaving only one
		# integration per group
		to_remove = []
		for group in same:
			
			# on each loop, add one alignment to to_remove until there is only one left
			while len(group) > 1:
				removed = False
				# preferentially remove any alignments that have hard clips
				for i in range(len(group)):
					if i in hard_clips:
						to_add = group.pop(i)
						to_remove.append(to_add)
						removed=True
						break
					
				
				# if no hard clips, just remove the last integration
				if not removed:
					to_add = group.pop(-1)
					to_remove.append(to_add)
				
		# now remove integrations
		assert len(to_remove) == len(set(to_remove))
		for i in reversed(sorted(to_remove)):
			self.pop(i)		
			
	def remove_query_nested(self):
		""" 
		Remove any alignments that are encompassed by other alignments.  That is,
		remove any alignments that cover the same part of the query, keeping the alignment 
		that covers the larger part of the query
		
		For example, considering two alignments, one 150M and 130M20S, retain only 150M.
		
		For two alignments, 50M100S and 40S110M, keep both because neither alignment
		completely encompasses the other
		
		If two alignments are equivalent (i.e. cover same part of read), remove the one
		with the higher edit distance
		
		Note that this changes the order of alignments in the pool
		"""
		
		if len(self) < 2:
			return
			
		assert all([read.query_name == self[0].query_name for read in self[1:]])
		
		# sort alignments in order of query alignment start and query alignment end
		self.sort(key = lambda elem: (elem.query_alignment_start, elem.query_alignment_end))
		
		# look for nested alignments
		i = 1
		while i  < len(self):

			# if we have the same start
			if self[i].query_alignment_start == self[i-1].query_alignment_start:
			
				# if alignments cover the same part of the read
				if self[i].query_alignment_end == self[i-1].query_alignment_end:

					# remove the one with the higher edit distance
					try:
						if self[i].get_tags('NM') > self[i-1].get_tags('NM'):
							self.pop(i)
							continue
						else:
							self.pop(i-1)
							continue
					# if there's no edit distance, just remove the ith one
					except KeyError:
						self.pop(i)
						continue
					
				# if stop of (i) is equal to or before
				elif self[i].query_alignment_end > self[i-1].query_alignment_end:
					self.pop(i-1)
					continue

			# if start of ith alignment is after start of (i-1)th alignment	(because sorted)
			else:
				# but end is before, then nested
				if self[i].query_alignment_end <= self[i-1].query_alignment_end:
					self.pop(i)
					continue
			
			i += 1
			
	def _convert_cigarstring_to_cigartuples(self, cigar):
		""" Convert a CIGAR string to cigar tuples """
		
		assert re.match("^(\d+[MIDNSHP=X])+$", cigar)
		
		tuples = []
		ops = { 'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}
		for op in re.findall("\d+[MIDNSHP=X]", cigar):
			tuples.append((ops[op[-1]], int(op[:-1])))
			
		return tuples
		
	def _convert_hard_clips_to_soft(self, cigartuples):
		""" Return the cigartuples with hard clips converted to soft clips """

		for i in range(len(cigartuples)):
			if cigartuples[i][0] == 5:
				cigartuples[i] = (4, cigartuples[i][1])
				
		return cigartuples
		
	def _create_new_segment_from_primary(self, primary, header, rname, pos, ori, cigar,
											nm, mapq=None):
		""" 
		Create a new pysam.AlignedSegment from an existing (primary) pysam.AlignnedSegment
		"""
		
		a = pysam.AlignedSegment(header=header)
		a.query_name = primary.query_name
		if mapq is not None:
			a.mapping_quality = mapq
		
		assert ori == '-' or ori == '+'
		
		# if primary alignment has different orientation to this one, need to reverse complement
		if primary.is_reverse == (ori == '-'):
			a.query_sequence = primary.query_sequence
			a.query_qualities = primary.query_qualities
		else:
			a.query_sequence = self._reverse_complement(primary.query_sequence)
			quals = array.array('B', primary.query_qualities)
			quals.reverse()
			a.query_qualities = quals
			
		flag = 0x100 + 0x800
		if primary.is_paired:
			flag += 0x1
		if primary.is_proper_pair:
			flag += 0x2
		if primary.mate_is_unmapped:
			flag += 0x8
		if ori == '-':
			flag += 0x10
		if primary.mate_is_reverse:
			flag += 0x20
		if primary.is_read1:
			flag += 0x40
		if primary.is_read2:
			flag += 0x80
		
		a.flag = flag

		a.reference_id = header.get_tid(rname)
		a.reference_start = pos
		a.cigar = self._convert_cigarstring_to_cigartuples(cigar)
		
		a.next_reference_id = primary.next_reference_id
		a.next_reference_start = primary.next_reference_start 
		a.template_length = primary.template_length
		
		if nm is not None:
			a.tags = [('NM', nm)]
		else:
			a.tags = []

		return a	
		
	def _reverse_complement(self, seq):
		""" Return reverse complement """
		nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N' ,
			  'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n'}
		return ''.join(reversed([nn[i] for i in seq]))			
			
class ChimericIntegration:
	""" 
	A class to store/calculate the properties of a simple chimeric integration 
	(ie mapped/clipped and clipped/mapped) 
	"""
	
	def __init__(self, host, virus, host_header, virus_header, 
					host_sec, virus_sec,  tol, map_thresh=20, primary=True):
		""" Given a host and virus read that are chimeric, calculate the properties of the integration """

		self.map_thresh = map_thresh
		self.tol = tol
		self.primary = primary
		
		assert self.tol >= 0
		
		self.hread = self._combine_short_CIGAR_elements(host, host_header)
		self.hread = self._simplify_CIGAR(self.hread, host_header)

		self.vread = self._combine_short_CIGAR_elements(virus, virus_header)
		self.vread = self._simplify_CIGAR(self.vread, virus_header)

		# simple chimeric read has two CIGAR elements
		assert self._is_chimeric()
		
		self.hsec = self._process_non_primary_alignments(self.hread, host_sec, host_header)
		self.hhead = host_header
		
		self.vsec = self._process_non_primary_alignments(self.vread, virus_sec, virus_header)
		self.vhead = virus_header
		
		self.chr = host.reference_name
		self.virus = virus.reference_name

	def get_ambig_coords(self):
		""" Get coordinates of ambiguous sequence relative to read """
		
		# if matched part is first
		if self.hread.cigartuples[0][0] == 0:
			return sorted([
				self.hread.query_alignment_end,
				self.vread.query_alignment_start
				])
		else:
			return sorted([
				self.hread.query_alignment_start,
				self.vread.query_alignment_end
				])

	def get_ambig_seq(self):
		""" Get sequence of ambiguous base from read """

		ambig = self.get_ambig_coords()
		return self.hread.query_sequence[ambig[0]:ambig[1]]

	def gap_or_overlap(self):
		""" Return 'gap' if there is a gap between host and viral alignments, 'overlap'
		if there's an overlap, and 'none' otherwise """
		
		ambig = self.num_ambig_bases()
		
		# if total mapped bases is equal to read length, there's no overlap or gap
		if ambig == 0:
			return 'none'
		
		# if there's more mapped bases in both alignments, there's an overlap
		elif ambig > 0:
			return 'overlap'
		
		# if less, there's a gap	
		else:
			return 'gap'
		
	def get_host_coords(self):
		""" Get coordinates of ambiguous bases in host """
		
		return self._get_ambig_coords(self.hread)

	def get_host_edit_dist(self):
		""" Get the edit distance for the host alignment """
		
		try:
			return self.hread.get_tag("NM")
		except KeyError:
			return None

	def get_host_mapq(self):
		""" Get mapping quality for host alignment """
		return self.hread.mapping_quality

	def get_host_seq(self):
		""" Get the part of the read aligned to the host (excluding ambiguous bases) """
		return self.hread.query_alignment_sequence

	def get_integration_orientation(self):
		""" 
		Get if the junction is host/virus (left side of integration), or virus/host (right side)
		"""
		
		# if host alignment is matched, clipped
		if self.hread.cigartuples[0][0] == 0:
			return 'hv'
		
		# if host alignment is clipped, matched
		else:
			return 'vh'
	
	def get_integration_type(self):
		""" Return type of integration (chimeric or discordant) """
	
		return 'chimeric'
				
	def get_properties(self):
		""" 
		Return a dict containing the properties of the integration 
		Properties are:
		Chr: host chromosome on which integration occurs
		IntStart, IntStop: coordinates of ambiguous bases (gap or overlap) in host chrom
		VirusRef: virus which is integrated
		VirusStart, VirusStop: coordinates of ambiguous bases (gap or overlap) in virus 
		NoAmbiguousBase: Number of bases involved in gap or overlap
		OverlapType: Type of junction - gap, overlap or none (neither)
		HostSeq: Sequence of host portion of the read (including overlap if relevant)
		VirusSeq: Sequence of the viral portion of the read (including overlap if relevant)
		AmbiguousSeq: Sequence of ambiguous bases
		HostEditDist: Edit distance for host alignment
		ViralEditDist: Edit distance for viral alignment
		TotalEditDist: Sum of HostEditDist and ViralEditDist, plus number of bases in
		overlap (if present)
		PossibleHostTranslocation: True if read can be accounted for only by alignments
		to host genome with lower edit distance than integration
		PossibleVectorRearrangement: True if read can be accounted for only by alignments
		to virus/vector genome with lower edit distance than integration
		HostAmbiguousLocation: True if primary alignment has one or more equivalent (same
		edit distance and CIGAR) alignments elsewhere in the host genome
		ViralAmbiguousLocation: True if primary alignment has one or more equivalent (same
		edit distance and CIGAR) alignments elsewhere in the virus genome, or a different
		provided viral reference
		AlternativeInts: A list of the alternative integrations
		Type: 'chimeric' or 'discordant'
		HostMapQ: Mapping quality for host alignment
		ViralMapQ: Mapping quality for viral alignment
		ReadID: query name, with /1 appended for integrations only involving R1 and /2
		appended for integrations only involving R2
		ReadSeq: query sequence (may be reverse complemented if the read is alignned in
		reverse orientation in host alignment)
		
		Note that coordinates are 0-based
		Note that host/virus/ambig/read sequences might be reverse-complemented self.ints
		compared to original read if host or virus alignment was in reverse orientation
		"""
		alt_ints = self._get_alt_ints()
		host_coords = self.get_host_coords()
		virus_coords = self.get_viral_coords()
		
		return {
			'Chr': self.chr,
			'IntStart': host_coords[0],
			'IntStop': host_coords[1],
			'VirusRef': self.virus,
			'VirusStart': virus_coords[0],
			'VirusStop': virus_coords[1],
			'NoAmbiguousBases': self.num_ambig_bases(absolute = True),
			'OverlapType': self.gap_or_overlap(),
			'Orientation': self.get_integration_orientation(),
			'VirusOrientation': self.get_viral_orientation(),
			'HostSeq': self.get_host_seq(),
			'VirusSeq': self.get_viral_seq(),
			'AmbiguousSeq': self.get_ambig_seq(),
			'HostEditDist': self.get_host_edit_dist(),
			'ViralEditDist': self.get_viral_edit_dist(),
			'TotalEditDist': self.get_total_edit_dist(),
			'PossibleHostTranslocation': self.is_possible_translocation(),
			'PossibleVectorRearrangement': self.is_possible_virus_rearrangement(),
			'HostAmbiguousLocation': self.is_host_ambig_loc(),
			'ViralAmbiguousLocation': self.is_viral_ambig_loc(),
			'AlternativeInts': ";".join([int.short_str() for int in alt_ints]),
			'Type': self.get_integration_type(),
			'HostMapQ': self.get_host_mapq(),
			'ViralMapQ': self.get_viral_mapq(),
			'ReadID': self.get_read_id(),
			'ReadSeq': self.hread.query_alignment_sequence
		}

	def get_read_id(self):
		"""  Return the query name, with /1 appended for integrations only involving R1 
		and /2 appended for integrations only involving R2
		"""
		if self.hread.is_read1:
			return self.hread.query_name + "/1"
		elif self.hread.is_read2:
			return self.hread.query_name + "/2"
		else:
			return self.hread.query_name
	
	def get_total_edit_dist(self):
		""" 
		Get total edit distance, which is equal sum of host edit distance, viral edit
		distance, and number of bases in gap (if present)
		"""
		
		total_edit_dist = self.get_host_edit_dist() + self.get_viral_edit_dist()
		if self.gap_or_overlap() == 'gap':
			return total_edit_dist + abs(self.num_ambig_bases())
		else:
			return total_edit_dist 
		
	def get_viral_coords(self):
		""" Get coordinates of ambiguous bases in virus """
		
		return self._get_ambig_coords(self.vread)		
		
	def get_viral_edit_dist(self):
		""" Get the edit distance for the viral alignment """
		
		try:
			return self.vread.get_tag("NM")
		except KeyError:
			return None
	
	def get_viral_mapq(self):
		""" Get mapping quality for viral alignment """
		return self.vread.mapping_quality
			
	def get_viral_orientation(self):
		""" Return + for virus inserted in forward orientation, - for reverse """
		
		if self.hread.is_reverse:
			if self.vread.is_reverse:
				return '+'
			else:
				return '-'
		else:
			if self.vread.is_reverse:
				return '-'
			else:
				return '+'
		
	def get_viral_seq(self):
		""" Get the part of the read aligned to the virus (excluding ambiguous bases)"""
		
		return self.vread.query_alignment_sequence
		
	def num_ambig_bases(self, absolute = False):
		""" 
		Get the number of ambiguous bases.  This is equal to the number of bases accounted 
		for in both alignments minus the read length.  
		
		if absolute is False:
			Positive => overlap
			Negative => gap
			Zero => neither
		if absolute is True, just return the number of ambiguous bases
		"""
		
		# make sure that read length is the same for host and viral reads
		assert self.hread.infer_read_length() == self.vread.infer_read_length()
		
		# get total number of bases mapped in both alignments
		total_mapped = self.hread.query_alignment_length + self.vread.query_alignment_length
		
		# ambig bases is total mapped - read length
		ambig = total_mapped - self.hread.infer_read_length()
		
		if absolute:
			return abs(ambig)
		else:
			return ambig
		
	def is_host_ambig_loc(self):
		""" 
		Check if an integration has ambiguous location in host genome (i.e. host part of
		read has multiple equivalent alignments to different parts of host genome)
		"""
		return self._is_ambiguous_location(self.hread, self.hsec)
	
	def is_possible_translocation(self):
		""" 
		A read is a possible host transloaction if more bases in the read can be 
		accounted for by alignments to the host than would be if the read was a host/virus
		chimera
		"""

		return self._is_possible_rearrangement(self.hread, self.hsec)
	
	def is_possible_virus_rearrangement(self):
		""" 
		A read is a possible viral rearrangement if more bases in the read can be 
		accounted for by alignments to the virus than would be if the read was a host/virus
		chimera
		"""

		return self._is_possible_rearrangement(self.vread, self.vsec)
		
	def is_viral_ambig_loc(self):
		"""
		Check if an integration has ambiguous location in viral genome (i.e. viral part of
		read has multiple equivalent alignments to different parts of viral genome)
		"""	
		return self._is_ambiguous_location(self.vread, self.vsec)
		
	def short_str(self):
		"""
		Return a shorter string representation of the integration
		chr:start-stop,ori,virus:start-stop,vori,n_ambig,overlap_type,nm
		"""
		chr = self.chr
		hstart, hstop = self.get_host_coords()
		ori = self.get_integration_orientation()
		virus = self.virus
		vstart, vstop = self.get_viral_coords()
		vori = self.get_viral_orientation()
		n_ambig = self.num_ambig_bases(absolute = True)
		overlap_type = self.gap_or_overlap()
		nm = self.get_total_edit_dist()
		
		props = (
			f"{chr}:{hstart}-{hstop}",
			ori,
			f"{virus}:{vstart}-{vstart}",
			vori,
			str(n_ambig),
			overlap_type,
			str(nm)
		)
		
		return ",".join(props)

	def _combine_short_CIGAR_elements(self, read, header):
		"""
		Sometimes, short insertions, deletions and soft-clipped regions can cause a read
		to be missed because it doesn't meet our rather strict criteria for a chimeric read
		
		For example, viral CIGAR 1S100M49S and host CIGAR 101S49M would be missed as a
		potential integration because the viral cigar is clipped on both ends,
		even though it really looks like an integration.  
		
		To address this issue, combine any non-mapped short CIGAR elements ()
		that are shorter in length than a tolerance value (self.tol)
		
		Similar in spirity to '_simplify_CIGAR', but takes a slightly different approach (
		for example, deals differently with reads that are soft-clipped on both ends)
		"""
		
		if read.cigartuples is None:
			return read
			
		if len(read.cigartuples) == 1:
			return read
		
		# first look for elements to simplify
		curr_non_mapped = 0
		curr = []
		last_mapped = None
		
		# key is matched region to combine into, and value
		# is a list of elements to combine
		to_delete = []
		
		for i in range(len(read.cigartuples)):
			# if mapped
			if read.cigartuples[i][0] == 0:
				# check if curr_non_mapped is more than zero and less than tol
				if curr_non_mapped > 0 and curr_non_mapped <= self.tol:
					to_delete.append((i, curr))
				curr_non_mapped = 0
				curr = []
				last_mapped = i
			# if not mapped
			else:
				curr.append(i)
				curr_non_mapped += read.cigartuples[i][1]
		
		# check last element
		if curr_non_mapped > 0 and curr_non_mapped <= self.tol:
			to_delete.append((last_mapped, curr))
		
		if len(to_delete) == 0:
			return read

		# next, combine short elements with nearest matched region
		tmp_cigartuples = list(read.cigartuples)
		nm_offset = 0
		
		for elem in reversed(to_delete):

			i_matched = elem[0]
			assert read.cigartuples[i_matched][0] == 0
			assert tmp_cigartuples[i_matched][0] == 0
			i_del = elem[1]
			
			# if cigar operation consumes query, need to add to matched region			
			matched = tmp_cigartuples[i_matched][1]
			for idx in i_del:
				# I (1) and S (4) operations consume query
				if read.cigartuples[idx][0] in {1, 4}:
					matched += read.cigartuples[idx][1]
				
				# if this is a soft-clip, we also need to add it to the edit distance offset
				if read.cigartuples[idx][0] == 4:
					nm_offset += read.cigartuples[idx][1]
			
			# update mapped length
			tmp_cigartuples[i_matched] = (0, matched)
			
			# remove elements
			for idx in i_del:
				tmp_cigartuples.pop(idx)
		
		# check that we don't now have two matched regions next to each other
		if len(tmp_cigartuples) > 1:
			for i in reversed(range(len(tmp_cigartuples)-1)):
				if tmp_cigartuples[i][0] == 0:
					if tmp_cigartuples[i+1][0] == 0:
						matched = tmp_cigartuples[i][1] + tmp_cigartuples[i+1][1]
						tmp_cigartuples.pop(i+1)
						tmp_cigartuples[i] = (0, matched)	
			
			
		# if we're removing elements that consume the query and are before the first
		# mapped region for a forward read or after the last mapped region for a reverse 
		# read then we need to also adjust the mapping position
		deleted_inds = [tup[1] for tup in to_delete]
		#https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-a-list-of-lists
		deleted_inds = [item for sublist in deleted_inds for item in sublist]
		mapped = [ind for ind in range(len(read.cigartuples)) if read.cigartuples[ind][0] == 0]

		if read.is_reverse:
			# check if there are any regions we deleted from after last mapped region
			offset_op = range(mapped[-1], len(read.cigartuples))
			offset_op = [i for i in offset_op if i in deleted_inds]
			
			if len(offset_op) > 0:
			
				# only need to consider those that consume query
				pos_offset = [read.cigartuples[ind] for ind in offset_op]
				pos_offset = [tup[1] for tup in pos_offset if tup[0] in (1, 4)]
				pos_offset = -sum(pos_offset)
				
			else:
				pos_offset = 0
		else:
			# check if there are any regions we deleted from before first mapped region
			offset_op = range(mapped[0])
			offset_op = [i for i in offset_op if i in deleted_inds]
		
			if len(offset_op) > 0:
				# only need to consider those that consume query
				pos_offset = [read.cigartuples[ind] for ind in offset_op]
				pos_offset = [tup[1] for tup in pos_offset if tup[0] in (1, 4)]
				pos_offset = -sum(pos_offset)
				
			else:
				pos_offset = 0


		return self._edit_alignment(read, header, tmp_cigartuples, pos_offset, nm_offset)
		
	def _edit_alignment(self, read, header, cigartuples, pos_offset, nm_offset):
		""" Edit an alignment by modifying the cigartuples of an existing read """
		a = pysam.AlignedSegment(header=header)
		a.query_name = read.query_name
		a.query_sequence = read.query_sequence
		a.flag = read.flag
		a.reference_id = read.reference_id
		a.reference_start = read.reference_start + pos_offset
		a.mapping_quality = read.mapping_quality
		a.cigar = cigartuples
		a.next_reference_id = read.next_reference_id
		a.next_reference_start = read.next_reference_start 
		a.template_length = read.template_length
		a.query_qualities = read.query_qualities
		
		# update edit distance in NM tag (but not MD - this is not reliable!!)
		tags = read.get_tags()
		updated_OA = False
		orig_strand = '-' if read.is_reverse else '+'
		OA_update = (read.query_name, str(read.reference_start), orig_strand, 
						read.cigarstring, str(read.mapping_quality), 
						str(read.get_tag('NM')))
		OA_update = ",".join(OA_update) + ";"
		for i in range(len(tags)):
			if tags[i][0] == 'NM':
				tags[i] = ('NM', tags[i][1] + nm_offset)
			if tags[i][0] == 'OA':
				tags[i] = ('OA', tags[i][1] + OA_update)
				updated_OA = True
		
		if not updated_OA:
			tags.append(('OA', OA_update))
			
		a.tags = tags
		
		return a
		
	def _get_alt_ints(self):
		"""
		If there is one or more secondary/supplementary alignment(s) that cover(s) the same 
		part of the read as the primary alignment, then we say that the location of the
		integration is ambiguous (since we don't know where the integration actually 
		originated from).  In the most complex case, we consider each viral alignment
		and each host alignment as constituting a possible alternative integration
		
		Return a list of ChimericIntegration objects with alternative integrations
		"""
		
		alt_ints = []
		
		host_alns = AlignmentPool(self.hsec)
		host_alns.append(self.hread)
		host_alns.remove_redundant_alignments()
		
		virus_alns = AlignmentPool(self.vsec)
		virus_alns.append(self.vread)
		virus_alns.remove_redundant_alignments()
			
		for hst in host_alns:
			for vrs in virus_alns:
				
				# don't do primary vs primary
				if hst == self.hread and vrs == self.vread:
					continue
					
				# if hst and vrs constitute a valid integration, add to list
				try:
					alt_int = ChimericIntegration(hst, vrs, self.hhead, self.vhead,
					AlignmentPool(), AlignmentPool(), self.tol, self.map_thresh, False)
					alt_ints.append(alt_int)
					
				except AssertionError:
					pass
			
		# if this integration is a primary one, but we found an alternative integration with
		# a lower edit distance, then the alternate one becomes primary
		primary_nm = self.get_total_edit_dist()
		if self.primary:
		
			for i in range(len(alt_ints)):
				
				alt_nm = alt_ints[i].get_total_edit_dist()
			
				if alt_nm < primary_nm:
			
					print('found better integration than primary one: swapping')
				
					# add secondary alignments to alt_int (not added during instatiation)
					host_alns.append(self.hread)
					host_alns.remove(alt_ints[i].hread)
				
					virus_alns.append(self.vread)
					virus_alns.remove(alt_ints[i].vread)
				
					alt_ints[i].hsec = host_alns
					alt_ints[i].vsec = virus_alns
				
					# swap self with alt_int
					tmp = self
					self = alt_ints[i]
					self.primary = True
					alt_ints[i] = tmp
					alt_ints[i].primary = False	

		return alt_ints
			
	def _get_ambig_coords(self, read):
		""" Return a tuple of (start, stop) of ambiguous bases in host or virus """
		# if aligned part of read is first
		if read.cigartuples[0][0] == 0:
			if self.gap_or_overlap() == 'gap':
				return (
					read.get_blocks()[0][1],
					read.get_blocks()[0][1] + abs(self.num_ambig_bases())
				)
			elif self.gap_or_overlap() == 'overlap':
				return (
					read.get_blocks()[0][1] - abs(self.num_ambig_bases()),
					read.get_blocks()[0][1]
				)
			else:
				return (
					read.get_blocks()[0][1],
					read.get_blocks()[0][1]
				)
		# if aligned part of read is second
		else:
			if self.gap_or_overlap() == 'gap':
				return (
					read.get_blocks()[0][0],
					read.get_blocks()[0][0] + abs(self.num_ambig_bases())
				)
			elif self.gap_or_overlap() == 'overlap':
				return (
					read.get_blocks()[0][0] - abs(self.num_ambig_bases()),
					read.get_blocks()[0][0]
				)
			else:
				return (
					read.get_blocks()[0][0],
					read.get_blocks()[0][0]
				)
				
	def _is_ambiguous_location(self, read, sec):
		"""  
		An alignment is ambiguous if there's a secondary or supplementary alignment
		that is equivalent to the primary one.
		An equivalent alignment is one with the same CIGAR and edit distance
		"""	

		if len(sec) == 0:
			return False
			
		try:
			primary_nm = read.get_tag('NM')
		except KeyError:
			primary_nm = None
	
		for sec_read in sec:
		
			# check for the same CIGAR
			same_cigar = sec_read.cigartuples == read.cigartuples

			# CIGAR must be the same
			if not same_cigar:
				continue
				
			# check for same edit distance
			try:
				sec_nm = sec_read.get_tag('NM')
			except KeyError:
				sec_nm = None
			
			if primary_nm is not None and sec_nm is not None:
				if primary_nm != sec_nm:
					continue 
				
			# check for different mapping position or strand
			same_pos = sec_read.reference_start == read.reference_start
			same_dir = sec_read.is_reverse == read.is_reverse
						
			# if position or strand is different, then location is ambiguous
			if (not same_pos) or (not same_dir):
				return True

		# if we didn't find any equivalent alignments, location not ambiguous
		return False
		
	def _is_possible_rearrangement(self, primary, sec):
		""" 
		Given a primary alignment and an AlignmentPool of secondary alignments, compare 
		edit distances for the case that the read is an integration, and the case that the
		read is a rearrangement
		
		Return True if rearrangement is more likely, otherwise False
		"""
		
		# if no secondary alignments, not possible rearrangement
		if len(sec) == 0:
			return False
			
		# make alignment pool with secondary and primary alignments
		all = AlignmentPool(sec)
		all.append(primary)
		
		rearrange_nm, all = all.get_rearrangement_nm()
		
		# if there's only one alignment in pool after removing redundant alignments,
		# can't be rearrangement
		if len(all) == 1:
			return False
			
		# compare integration edit distance with rearrangement edit distance
		int_nm = self.get_total_edit_dist()
		
		if rearrange_nm <= int_nm:
			return True
		else:
			return False
					
	def _is_chimeric(self):
		"""
		Return True if two alignments are chimeric (occupy complementary parts of the read)
		Mapped parts of the read must be at least self.map_thresh bp long
		
		"""
		
		# if not aligned, can't be chimeric
		if self.vread.is_unmapped or self.hread.is_unmapped:
			return False
			
		# viral and host alignment must have only mapped and clipped regions
		if not (len(self.vread.cigartuples) == 2 and len(self.hread.cigartuples) == 2):
			return False
		
		# https://pysam.readthedocs.io/en/latest/api.html#api
		# 0 = matched
		# 4 = soft-clipped
		cigarops1 = [i[0] for i in self.hread.cigartuples]
		cigarops2 = [i[0] for i in self.vread.cigartuples]
		
		# must have one matched and one soft-clipped part
		if set(cigarops1) != {0, 4}:
			return False
		if set(cigarops2) != {0, 4}:
			return False
			
		# mapped regions must be at least 20 bp
		if self.hread.query_alignment_length < self.map_thresh:
			return False
		if self.vread.query_alignment_length < self.map_thresh:
			return False
			
		# if both forward or both reverse
		if self.hread.is_reverse == self.vread.is_reverse:
			# but cigar operations are the same, then no dice
			if cigarops1 == cigarops2:
				return False
		# if aligned in opposite orientations
		else:
			if cigarops1 != cigarops2:
				return False	

		return True		
		
	def _process_non_primary_alignments(self, primary, not_primary, header):
		""" 
		Depending on aligner, secondary and supplementary alignments might be present in
		tags from the primary alignment, or as separate alignments, or both.  
		
		We might also want to simplify these alignments using self._combine_short_CIGAR_elements()
		and self._simplify_CIGAR()
		"""
		
		not_primary = self._process_SA_and_XA(primary, header)
		
		# combine any short elements in secondary alignments
		for i in range(len(not_primary)):
			not_primary[i] = self._combine_short_CIGAR_elements(not_primary[i], header)
		
		if len(not_primary) < 2:
			return not_primary
			
		# remove any redundant alignments from pool
		# redundant alignments have same CIGAR, orientation (+/-) and pos		
		not_primary.remove_redundant_alignments()
		
		return not_primary
	
	def _process_SA(self, SA, primary, header):
		""" 
		Create pysam.AlignedSegment objects for each XA alignment 
		
		SA has format (rname,pos,strand,CIGAR,mapQ,NM;)+
		"""
		
		SA_segs = AlignmentPool()
		
		SA = SA.split(";")
		SA = [i for i in SA if i != '']
		
		for aln in SA:
			SA_segs.add_read_from_SA(aln, primary, header)
		
		return SA_segs
				
	def _process_SA_and_XA(self, primary, header):
		""" Create AlignmentPool objects with alignments from XA and SA """
		
		non_primary = AlignmentPool()
		
		# get supplementary alignments from SA tag
		try:
			SA = primary.get_tag('SA')
			non_primary += self._process_SA(SA, primary, header)
		except KeyError:
			pass
		
		# get supplementary alignments from XA tag
		try:
			XA = primary.get_tag('XA')
			non_primary += self._process_XA(XA, primary, header)
		except KeyError:
			pass
			
		return non_primary		
	
	def _process_XA(self, XA, primary, header):
		""" 
		Create pysam.AlignedSegment objects for each XA alignment 
		
		XA is bwa-specific and has format (chr,pos,CIGAR,NM;)*
		"""
		
		XA = XA.split(";")
		XA = [i for i in XA if i != '']
		
		XA_segs = AlignmentPool()
		for aln in XA:
			XA_segs.add_read_from_XA(aln, primary, header)			
		
		
		return XA_segs		
		
	def _simplify_CIGAR(self, read, header):
		""" 
		Simplify CIGAR by combining all matched regions into one, leaving any
		soft-clipped regions at either end intact
		
		When combining, consider operations that do and don't consume the query
		when deciding how many bases should be added to the matched region (I - consumes query;
		DNP - don't consume query)
		
		For example: 
		10S40M30I40M => 10S110M
		10S40M20D40M => 10S80M
		10S40M30I20D40M => 10S110M
		
		This is similar to _combine_short_CIGAR_elements, but can combine larger operations
		that break up a matched region
		"""
		
		if read.cigartuples is None:
			return read
			
		matched = [i for i in range(len(read.cigartuples)) if read.cigartuples[i][0] == 0]
		if len(matched) == 1:
			return read
		
		# figure out how long the new matched region should be
		n_match = 0
		tmp_cigartuples = read.cigartuples
		for i in reversed(range(matched[0], matched[-1]+1)):
			if read.cigartuples[i][0] in (0, 1):
				n_match += read.cigartuples[i][1]
			tmp_cigartuples.pop(i)
			
		# insert just one matched region
		tmp_cigartuples.insert(matched[0], (0, n_match))
				
		return self._edit_alignment(read, header, tmp_cigartuples, 0, 0)	
		
	def __str__(self):
		
		props = self.get_properties()
		
		return "\t".join([str(props[i]) for i in default_header])

class DiscordantIntegration(ChimericIntegration):
	""" 
	A class to store/calculate the properties of a discordant integration
	(i.e. read 1 mapped to one reference, read2 mapped to a different reference)
	
	This class inherits methods from ChimericIntegration, but re-defines methods that
	need to be implemented differently for a discordant pair
	"""
	
	def __init__(self, host_r1, host_r2, virus_r1, virus_r2, tol,
					host_header, virus_header, 
					host_sec1, host_sec2, virus_sec1, virus_sec2,
					map_thresh=20, tlen=0, primary=True):

		self.map_thresh = map_thresh
		self.tlen = tlen
		self.tol = tol
		self.host_header = host_header
		self.virus_header = virus_header
		self.primary = primary
		
		self.hread1 = self._combine_short_CIGAR_elements(host_r1, host_header)
		self.hread2 = self._combine_short_CIGAR_elements(host_r2, host_header)
		self.vread1 = self._combine_short_CIGAR_elements(virus_r1, virus_header)
		self.vread2 = self._combine_short_CIGAR_elements(virus_r2, virus_header)
		
		self.hsec1 = host_sec1
		self.hsec2 = host_sec2
		self.vsec1 = virus_sec1
		self.vsec2 = virus_sec2	
					
		assert tol >= 0
		assert tlen >= 0
		
		# note that self._is_discordant() also assigns self.hread and self.vread
		assert self._is_discordant(host_sec1, host_sec2, virus_sec1, virus_sec2)
		
		self.get_host_chr()
		self.get_viral_ref()
		
		
	def gap_or_overlap(self):
		""" For a discordant pair, OverlapType is always discordant """
		return 'discordant'
	
	def get_ambig_seq(self):
		""" No ambiguous bases for a discordant pair"""

		return None

	def get_coords(self, read, header):
		"""
		Get the likely coordinates of integration.  This is the end of the aliged read
		that faces the other read, plus the mean insert size
		"""


		# mean insert length is mean template length minus sum of R1 and R2 lengths
		insert = self.tlen - self.vread.query_length - self.hread.query_length
		# but can't be less than zero
		if insert < 0:
			insert = 0
		
		# if host read is R1, then we want the coordinate of the host base at the
		# side of the read that faces R2
		if read.is_read1:
			if read.is_reverse:
				stop = read.get_blocks()[0][0]
				start = stop - insert
			else:
				start = read.get_blocks()[-1][-1]
				stop = start + insert
		else:
			if read.is_reverse:
				start = read.get_blocks()[-1][-1]
				stop = start + insert
			else:
				stop = read.get_blocks()[0][0]
				start = stop - insert
			
		# check the start isn't less than zero
		if start < 0:
			start = 0

		# check that the stop isn't more than the length of the reference
		reference_length = header.get_reference_length(read.reference_name)
		if stop > reference_length:
			stop = reference_length
		
		return (start, stop)

	def get_host_coords(self):
		""" Get likely coordinates of integration in host chromosome """
		
		return self.get_coords(self.hread, self.host_header)

	def get_host_chr(self):
		""" Get reference name for host alignment """
		
		h1_map, h2_map = self._is_mapped(self.hread1), self._is_mapped(self.hread2)
		
		assert h1_map != h2_map
		
		if h1_map:
			self.chr = self.hread1.reference_name
		else:
			self.chr = self.hread2.reference_name			
			
	def get_integration_orientation(self):
		""" 
		Orientation refers to side of the integration that we see in the read: 
		host/virus (i.e. left side of the integration) or virus/host (i.e. right side)
		"""

		if self.hread.is_reverse:
			return 'vh'
		else:
			return 'hv'
	
	def get_integration_type(self):
		""" Get type of integration - in this case discordant """
	
		return 'discordant'
	
	def num_ambig_bases(self, absolute=True):
		""" 
		For a discordant pair, the number of ambiguous bases is not defined.
		Need to keep the absolute argument for consistency with ChimericIntegration
		"""
		
		return None
	
	def get_read_id(self):
		""" Return the query name / read ID """
		return self.hread.query_name
		
	def get_total_edit_dist(self):
		""" 
		The total edit distance reflects edits from a 'perfect' integration.  In the context
		of a chimeric read, this is one in which the host and viral alignments have an 
		edit distance of 0, and there is no gap between the host and viral alignments.
		However, in the context of a discordant integration, the situation is a bit more
		complex.
		
		Ideally, we should consider soft-clipped bases when calculating an edit distance.
		
		For example, if we have the following (o = mapped, x = soft-clipped)
		              read 1                 read 2
		host:  xxxxxxxxxxxxxxxxxxxxx ooooooooooooooooooooo
		virus: xxxoooooooooooooooxxx xxxxxxxxxxxxxxxxxxxxx
		       ^^^               ^^^
		Should we consider the bases marked with ^ in the edit distance?
		
		Soft-clipped bases that are on the far side of the read from the other read
		(e.g. those marked ^ on the left in the example above) are unexpected, because in the simplest
		case we expect one read to be entirely mapped and the other entirely unmapped.  So
		perhaps we should add the number of these bases to the edit distance.
		
		Soft-clipped bases that are on the side of the read that faces the other read (e.g.
		those marked ^ on the right in the example above) may or may not be unexpected.
		In the example above, they are unmapped in both alignments, so perhaps they should
		be added to the edit distance.  However, also consider a situation like this:
		
		              read 1                 read 2
		host:  xxxxxxxxxxxxxxxxxxooo ooooooooooooooooooooo
		virus: ooooooooooooooooooxxx xxxxxxxxxxxxxxxxxxxxx
                                 ^^^
        In this case, these soft-clipped bases are probably not unexpected - they indicate 
        that the host/virus junction lies just before the end of read 1.  In practice, 
        we are not likely to actually see evidence of arrangement though, because there is 
        a minimum length of alignments (e.g. the seed length of the aligner).  For example,
        if the junction is 5bp from the right of read 1, all we will see is that the
        viral alignment is soft-clipped (5bp), but since 5bp is smaller than our seed length
        we won't see the 5bp alignment in the host alignment for read 1.  So we might 
        conclude that we should add 5 to the edit distance for the 5 soft-clipped bases
        on the right end of the viral alignment for read 1, but we'd probably be wrong.
        
        We therefore don't have enough information to calculate an edit distance analogous
        to the chimeric case for a discordant integration (most of the time). We could
        output the sum of the host and viral alignment edit distances, but I think this
        is misleading because it doesn't (and can't) capture the whole picture.  So instead,
        just output None
		"""
		return None
	
	def get_viral_coords(self):
		""" Get likely coordinates of integration in viral reference """		

		return self.get_coords(self.vread, self.virus_header)
		
	def get_viral_orientation(self):
		""" 
		Figure out if the virus is integrated in forward (+) or reverse (-) orientation
		"""

		if self.get_integration_orientation() == 'vh':
			if self.vread.is_reverse:
				return '-'
			else:
				return '+'	
		else:
			if self.vread.is_reverse:
				return '+'
			else:
				return '-'						
			
	def get_viral_ref(self):
		""" Get reference name for viral alignment """
		
		v1_map, v2_map = self._is_mapped(self.vread1), self._is_mapped(self.vread2)
		
		assert v1_map != v2_map
		
		if v1_map:
			self.virus = self.vread1.reference_name
		else:
			self.virus = self.vread2.reference_name
					
	def _is_discordant(self, host_sec1, host_sec2, virus_sec1, virus_sec2):
		""" 
		Check if a read-pair is discordant integration (one read mapped to host, one
		read mapped to virus)
		"""
		h1_map = self._is_mapped(self.hread1)
		h2_map = self._is_mapped(self.hread2)
		v1_map = self._is_mapped(self.vread1)
		v2_map = self._is_mapped(self.vread2)
		
		# pair from each alignment must have one mapped and one unmapped
		if h1_map == h2_map:
			return False
		if v1_map == v2_map:
			return False
			
		# must have complementary reads mapped (i.e. one that's mapped in one alignment
		# must be unmapped in the other)
		if h1_map != v2_map:
			return False
		if h2_map != v1_map:
			return False
			
		# assign self.hread and self.vread to help with some methods inherited from ChimericIntegration
		if h1_map:
			self.hread = self.hread1
			self.vread = self.vread2
			self.hsec = self._process_non_primary_alignments(self.hread, host_sec1, self.host_header)
			self.vsec = self._process_non_primary_alignments(self.vread, virus_sec2, self.virus_header)
		else:
			self.hread = self.hread2
			self.vread = self.vread1
			self.hsec = self._process_non_primary_alignments(self.hread, host_sec2, self.host_header)
			self.vsec = self._process_non_primary_alignments(self.vread, virus_sec1, self.virus_header)

		return True
		
	def _is_mapped(self, read):
		""" 
		For a discordant pair, a read is considered mapped if it has at most 20 unmapped 
		bases.  
		"""

		# if read is completely unmapped
		if read.is_unmapped:
			return False
			
		# if mapped pat of the read is less than threshold, is unmapped
		if read.query_alignment_length < self.map_thresh:
			return False
		
		# if non-mapped part of read is more than 20bp (or whatever threshold is), is unmapped
		if read.query_length - read.query_alignment_length > self.map_thresh:
			return False
					
		return True		
		
	def _is_possible_rearrangement(self, hread, vread):
		""" 
		For a DiscordantIntegration, rearrangements are unlikely since they would require 
		both reads to be mapped to the same reference, but a DiscordantIntegration is one
		in which one read in the pair is mapped and the other unmapped (when considering a
		single reference).  
		
		The only way it would be possible is if there are several 
		alignments which span the 'unmapped' read, each of which alone doesn't cover 
		enough of the read for it to be considered mapped, but when considered together
		they span most or all of the read.  Since we enforce a maximum number of mapped
		bases for the unmapped bases (by default 20), then all alignments would have to be
		shorter than this maximum.  So for a 150bp read, we would need at least 8 alignments
		of length 20bp to cover the whole read.  However, it's extremely unlikely that 8
		short fragments of the vector were recombined to create the unmapped read.  
		
		Therefore, assume that rearrangments are exceedingly unlikely for 
		DiscordantIntegrations and just return False.  
		"""	
		
		return False
		
	def _get_alt_ints(self):
		"""
		If there is one or more secondary/supplementary alignment(s) that cover(s) the same 
		part of the read as the primary alignment, then we say that the location of the
		integration is ambiguous (since we don't know where the integration actually 
		originated from).  In the most complex case, we consider each viral alignment
		and each host alignment as constituting a possible alternative integration.
		
		In particular, in the case of a DiscordantIntegration, all the properties are
		determined by the mapped read (not the unmapped read), so only consider
		the secondary and supplementary alignments for the mapped read for each reference 
		
		Return a list of DiscordantIntegration objects with alternative integrations
		"""

		alt_ints = []
		
		hread1_map = (self.hread1 == self.hread)
		
		if hread1_map:
			haln = AlignmentPool(self.hsec1)
			haln.append(self.hread1)
			haln.remove_redundant_alignments()
			
			valn = AlignmentPool(self.vsec2)
			valn.append(self.vread2)
			valn.remove_redundant_alignments()
			
			h2 = self.hread2
			v1 = self.vread1
		
		else:
			haln = AlignmentPool(self.hsec2)
			haln.append(self.hread2)
			haln.remove_redundant_alignments()
		
			valn = AlignmentPool(self.vsec1)
			valn.append(self.vread1)
			valn.remove_redundant_alignments()
			
			h1 = self.hread1
			v2 = self.vread2

		if len(haln) == 1 and len(valn) == 1:
			return alt_ints
			
		for host, virus in itertools.product(haln, valn):
				
				# don't do primary vs primary
				if host == self.hread and virus == self.vread:
					continue
					
				if hread1_map:
					h1, v2 = host, virus
				else:
					h2, v1 = host, virus

				# if host and virus reads constitute a valid integration, add to list
				try:
					alt_int = DiscordantIntegration(h1, h2, v1, v2, self.tol, 
					self.host_header, self.virus_header, AlignmentPool(), AlignmentPool(),  
					AlignmentPool(), AlignmentPool(), self.tol, self.map_thresh, False)
					alt_ints.append(alt_int)
					
				except AssertionError:
					pass

		return alt_ints
		
if __name__ == "__main__":
	main()