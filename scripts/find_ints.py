#!/usr/bin/env python3

import argparse
import pysam
import csv
import re
import array
import itertools
import copy
import math
import pdb

# note python 3.7+ is required

default_header = ['Chr', 'IntStart', 'IntStop', 'VirusRef', 'VirusStart', 'VirusStop', 
		'NoAmbiguousBases', 'OverlapType', 'Orientation', 'VirusOrientation', 
		'HostSeq', 'VirusSeq', 'AmbiguousSeq', 'HostEditDist', 'ViralEditDist', 
		'TotalEditDist', 'PossibleHostTranslocation', 'PossibleVectorRearrangement', 
		'HostAmbiguousLocation', 'ViralAmbiguousLocation', 'Type', 
		'HostMapQ', 'ViralMapQ', 'ReadID', 'AltLocs', 'ReadSeq']

def main():
	
	parser = argparse.ArgumentParser(description='Identify integrations by comparing a host and viral alignment')
	parser.add_argument('--host', help="host alignment bam (query-sorted)", required=True)
	parser.add_argument('--virus', help="virus alignment bam (query-sorted)", required=True)
	parser.add_argument('--integrations', help='file to output integrations to', default='integrations.tsv')
	parser.add_argument('--map-thresh', help="threshold for length of mapped region", default=20, type=int)
	parser.add_argument('--mean-template-length', help="mean template length for library", default=0, type=float)
	parser.add_argument('--tolerance', help="tolerance for short CIGAR operations (MIS ops with a combined length equal to or shorter than this will be combined into the nearest mapped region)", type=int, default=3)	
	parser.add_argument('--verbose', '-v', help="print extra messages?", action='store_true')
	parser.add_argument('--debug', '-d', help="print lots of extra messages?", action='store_true')
	parser.add_argument('--nm-diff', help="output all alternate integration sites with an edit distance at most this much more than the primary integration", type=int, default=None)
	parser.add_argument('--nm-pc', help="output all alternate integration sites with an edit distance at most this percentage more than the primary integration", type=float, default=None)
	args = parser.parse_args()

	with AlignmentFilePair(args.host, args.virus, args.integrations, args.map_thresh, args.mean_template_length,
							args.tolerance, args.nm_diff, args.nm_pc, args.verbose, args.debug) as pair:		
		
		pair.find_integrations()

		
class AlignmentFilePair:
	""" A class to hold host and viral alignments, and look for integrations """
	
	def __init__(self, host, virus, outfile, map_thresh, tlen, tol = 3, nm_diff=None, nm_pc=None, verbose = False, debug=False):
		""" Open host and viral alignments, and check for query sorting """
		
		self.host = AlignmentFile(host, verbose)
		self.virus = AlignmentFile(virus, verbose)	
		
		
		self.href_lens = self.host.get_reference_lengths()
		self.vref_lens = self.virus.get_reference_lengths()
		self.map_thresh = map_thresh
		self.verbose = verbose
		self.debug = debug
		self.tlen = int(tlen) # tlen needs to be int so that output coords are ints
		self.tol = tol
		self.nm_diff = nm_diff
		self.nm_pc = nm_pc
		
		# open outfile and write a header
		self.outfile = outfile
		self.out = open(outfile, 'w')
		self.header = default_header
		self.out.write("\t".join(self.header) + "\n")
		
		if self.tlen == 0:
			print(f"warning: insert size is zero - the position for discordant integration sites will be at the end of the mapped read")
		
		if self.tlen < 0:
			raise ValueError("Mean template length must be a positive integer")
		
		if self.tol < 0:
			raise ValueError("Tolerance must be a positive integer")
		
		if self.verbose:
			print(f"using host bam {host} and virus bam {virus}")
			

	def find_integrations(self):
		""" Check for integrations """
		
		counter = 0
		progress = 100000
		ints = 0
		
		# iterate over reads that have alignments to both host and virus
		for _ in self._get_aligns():
			
			if self.debug:
				print(f"{self.__repr__()}: find_integrations {self.curr_host[0].query_name}")
		
			counter += 1
			if counter % progress == 0:
				print(f"checked {counter} reads")
			
			if self.curr_host is None or self.curr_virus is None:
				break
				
			if self.verbose:
				print(f"checking {self.curr_host[0].query_name}")
				
			# collect read 1 alignments
			host_r1_primary  = self.curr_host.get_primary_r1(self.tol)
			host_r1_not_primary = self.curr_host.get_not_primary().get_read1()
			
			virus_r1_primary =  self.curr_virus.get_primary_r1(self.tol)
			virus_r1_not_primary = self.curr_virus.get_not_primary().get_read1()

			# check if primary read 1 alignments appear to be chimeric	
			if host_r1_primary is not None and virus_r1_primary is not None:
				if self._is_chimeric(host_r1_primary, virus_r1_primary,
									host_r1_not_primary, virus_r1_not_primary):
					ints += 1
												
				# or are a full integration
				if self._is_full(host_r1_primary, virus_r1_primary,
									host_r1_not_primary, virus_r1_not_primary):
					ints += 1

																				
			# collect read2 alignments
			host_r2_primary = self.curr_host.get_primary_r2(self.tol)
			host_r2_not_primary = self.curr_host.get_not_primary().get_read2()
			
			virus_r2_primary = self.curr_virus.get_primary_r2(self.tol)
			virus_r2_not_primary = self.curr_virus.get_not_primary().get_read2()			

			
			# check if primary read 1 alignments appear to be chimeric	
			if host_r2_primary is not None and virus_r2_primary is not None:
				if self._is_chimeric(host_r2_primary, virus_r2_primary,
									host_r2_not_primary, virus_r2_not_primary):
					ints += 1
		
				if self._is_full(host_r2_primary, virus_r2_primary,
									host_r2_not_primary, virus_r2_not_primary):
					ints += 1

										
			# if there are alignments that are neither read 1 nor read 2, check if chimeric		
			host_not12_primary = self.curr_host.get_primary_not_r1_or_r2(self.tol)
			host_not12_not_primary = self.curr_host.get_not_primary().get_not_read1_or_read2()
			
			virus_not12_primary = self.curr_virus.get_primary_not_r1_or_r2(self.tol)
			virus_not12_not_primary = self.curr_virus.get_not_primary().get_not_read1_or_read2()
			
			if host_not12_primary is not None and virus_not12_primary is not None:
				if self._is_chimeric(host_not12_primary, virus_not12_primary,
									host_not12_not_primary, virus_not12_not_primary):
					ints += 1
					
				if self._is_full(host_not12_primary, virus_not12_primary,
									host_not12_not_primary, virus_not12_not_primary):
					ints += 1
					
			
			# check for a discordant pair
			if all([i is not None for i in [host_r1_primary, host_r2_primary, 
											virus_r1_primary, virus_r2_primary]]):
				# check if integration
				if self._is_discordant(host_r1_primary, host_r2_primary, 
										virus_r1_primary, virus_r2_primary,
										host_r1_not_primary, host_r2_not_primary,
										virus_r1_not_primary, virus_r2_not_primary):
					ints += 1		
								
		print(f"found {ints} integrations")
		print(f"saved output to {self.outfile}")				
	
	def _write_integration(self, integration):
		"""
		Write an integration to file
		"""
		self.out.write(str(integration) + "\n") 
				
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
		if self.debug:
			print("\tchecking for chimeric integration")	
			
		try:
			integration = ChimericIntegration(hread, vread, map_thresh = self.map_thresh,
												tol = self.tol, nm_diff=self.nm_diff,
												nm_pc=self.nm_pc)	
		except AssertionError:
			return False

		if self.debug:
			print("\t setting secondary alignments")				
		integration.set_sec_alns(hsec, vsec)
		
		if self.debug:
			print("\t writing integration")		
		self._write_integration(integration)
		
		return True
			
	def _is_discordant(self, hread1, hread2, vread1, vread2, 
						hsec1, hsec2, vsec1, vsec2):
		"""
		Checks if a pair looks like an discordant integration
		"""
		if self.debug:
			print("\tchecking for discordant integration")				
		try:
			integration = DiscordantIntegration(hread1, hread2, vread1, vread2, 
										map_thresh = self.map_thresh, tlen = self.tlen,
										tol = self.tol, nm_diff=self.nm_diff,
										nm_pc=self.nm_pc)

		except AssertionError:
			return False

		if self.debug:
			print("\t setting secondary alignments")			
		integration.set_sec_alns(hsec1, hsec2, vsec1, vsec2)

		if self.debug:
			print("\t writing integration")			
		self._write_integration(integration)
		
		return True
			
	def _is_full(self, hread, vread, hsec, vsec):
		""" 
		Takes two alignments for the same read, one from the host and one from the virus, 
		and decides if the read appears to contain a 'full integration''.  A 'full 
		integration' means that one part of the read is accounted for by the host alignment, 
		the middle is accounted for by a viral alignment, and the other end is accounted 
		for by a host alignment.  For example 30M80I30M for host and 30S80M30S for virus.
		There might be some overlap or gap between the host and viral alignments, as
		for a ChimericIntegration
		"""					
		if self.debug:
			print("\tchecking for full integration")	
		try:
			integration = FullIntegration(hread, vread, map_thresh = self.map_thresh,
											tol = self.tol, nm_diff=self.nm_diff,
												nm_pc=self.nm_pc)

		except AssertionError:
			return False
		
		if self.debug:
			print("\t setting secondary alignments")
		integration.set_sec_alns(hsec, vsec)
			
		if self.debug:
			print("\t writing integration")
		self._write_integration(integration)
		
		return True

	def _get_aligns(self):
		""" A generator to get alignments for the same read from both host and virus"""

		# collect alignments from file
		self.curr_host = self.host.collect_next_read_alignments()
		self.curr_virus = self.virus.collect_next_read_alignments()
		
		prev_host = self.curr_host[0].query_name
		prev_virus = self.curr_virus[0].query_name

		# continue until we've either found the same read or reached the end of one of the files
		while True:
		
			if self.curr_host is None or self.curr_virus is None:
				return
				
			if not self._check_sort_order(self.curr_host, prev_host):
				raise ValueError(f"host bam must be sorted with 'samtools sort -n' (found out of order read {self.curr_host[0].query_name})")
			if not self._check_sort_order(self.curr_virus, prev_virus):			
				raise ValueError("host bam must be sorted with 'samtools sort -n' (found out of order read {self.curr_virus[0].query_name})")
				
			# check if the query names are the same for host and viral alignments
			if self.curr_host[0].query_name == self.curr_virus[0].query_name:
				
				yield
				self.curr_host = self.host.collect_next_read_alignments()
				self.curr_virus = self.virus.collect_next_read_alignments()

			# host query_name less than viral query_name
			elif self._str_gt(self.curr_host[0].query_name, self.curr_virus[0].query_name):
				if self.verbose:
					print(f"warning: host bam has no alignment for read {self.curr_virus[0].query_name} in virus bam, skipping")
				self.curr_virus = self.virus.collect_next_read_alignments()

			# viral query_name less than host query_name
#			elif self.curr_virus[0].query_name > self.curr_host[0].query_name:
			else:
				if self.verbose:
					print(f"warning: virus bam has no alignments for read {self.curr_host[0].query_name} in host bam, skipping")
				self.curr_host = self.host.collect_next_read_alignments()
				
			# take note of query name to check sort order
			try:
				prev_host = self.curr_host[0].query_name
			except TypeError:
				prev_host = None
			try:
				prev_virus = self.curr_virus[0].query_name
			except TypeError:
				prev_host = None
				
	def close(self):
		""" Close both alignments """

		self.host.close()
		self.virus.close()	
		self.out.close()
		
	def __enter__(self):
		return self
	
	def __exit__(self, exception, error, traceback):
		self.close()
		
	def _str_gt(self, str1, str2):
		"""
		Compare strings in the same way that samtools sorts: a 'natural' sort
		Break strings into chunks consisting of all letters and all numbers and compare these
		if the chunk contains non-digit characters sort lexiographically, otherwise sort 
		numerically
		"""
		
		#https://blog.codinghorror.com/sorting-for-humans-natural-sort-order/
		convert = lambda text: int(text) if text.isdigit() else text
		chunk = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
		
		for c1, c2 in zip(chunk(str1), chunk(str2)):
			if c1 == c2:
				continue
			if self._is_int(c1) and self._is_int(c2):
				return c1 > c2
			else:
				return str(c1) > str(c2)
		 
		return False
		
	def _is_int(self, x):
		try:
			int(x)
			return True
		except:
			return False
			
	def _check_sort_order(self, curr_lst, prev):
		"""
		Return True if current queryname is greater than or equal to previous queryname
		"""
		if curr_lst is None or prev is None:
			return True
		
		if curr_lst[0].query_name == prev:
			return True
		
		if self._str_gt(curr_lst[0].query_name, prev):
			return True
		
		return False

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
			
		query_name = self.curr.query_name
		
		if self.end:
			return None
		
		# keep getting elements until we either reach a different read or the end of the list
		while (self.curr.query_name == query_name) and (self.curr.query_name is not None):
			# includes R1, R2, primary, secondary, supplementary alignments
			alns.append(self.curr)		
			self._get_next_align()
			if self.end:
				break
				
		# check for reads that don't have a read sequence
		for read in alns:
			# if reads don't have a query sequence
			if read.query_length == 0:
				# try to get read sequence from another read
				for read2 in alns:
					if read == read2:
						continue
					
					if not (read.is_read1 == read2.is_read1 and read.is_read2 == read2.is_read2):
						continue
					if read2.query_length == 0:
						continue
					if read.is_reverse == read2.is_reverse: 
						read.query_sequence = str(read2.query_sequence)
#						read.query_qualities = copy.copy(read2.query_qualities)
					else:
						seq = _reverse_complement(str(read2.query_sequence))
						read.query_sequence = seq
#						read.query_qualities = copy.copy(read2.query_qualities).reverse()
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
				self.curr.query_name = None
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
	
	Should only contain alignments from the same read-pair (i.e. query_name should be the same for all reads in pool)
	"""
	
	def add_read_from_XA(self, XA, primary):
		""" Add a read based on an XA string (chr,pos,CIGAR,NM) """

		
		# get info from XA string
		XA = XA.split(",")
		assert len(XA) == 4
		
		chr = XA[0]
		ori = XA[1][0]
		pos = int(XA[1][1:]) - 1
		cigar = XA[2]
		nm = int(XA[3])
		
		# use this info with primary alignment to create new pysam.AlignedSegment 
		new_seg = self._create_new_segment_from_primary(primary, ref_name=chr, pos=pos, 
														ori=ori,cigar=cigar, nm=nm)
		
		# we copied all flags from primary, so need to update this												
		new_seg.is_supplementary = True
		# remove MD tag since it only applies to primary alignment
	
		new_seg = self._remove_MD_tag(new_seg)
		
		self.append(new_seg)
	
	def add_read_from_SA(self, SA, primary):
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
		
		new_seg = self._create_new_segment_from_primary(primary, ref_name=chr, pos=pos, 
														ori=ori, cigar=cigar, nm=nm, 
														mapq=mapq)
													
		new_seg.is_secondary = True

		new_seg = self._remove_MD_tag(new_seg)
		
		self.append(new_seg)
		
	def combine_short_CIGAR_elements(self, tol):
		"""
		combine short cigar elements in all reads in pool
		"""
		
		for i in range(len(self)):
			self[i] = self._combine_short_CIGAR_elements(tol, self[i])
	
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
		
	def get_primary_r1(self, tol):
		""" Return the primary read 1 alignment as a pysam.AlignedSegment object, or None if there isn't one """
		
		primary_r1 = self.get_primary().get_read1()
		if len(primary_r1) == 0:
			return None
		elif len(primary_r1) == 1:
			return self._combine_short_CIGAR_elements(tol, primary_r1[0])
		else:
			raise ValueError(f"found more than one primary alignment for {self[0].query_name} read 1")	
		
	def get_primary_r2(self, tol):
		""" Return the primary read 2 alignment as a pysam.AlignedSegment object, or None if there isn't one """
		
		primary_r2 = self.get_primary().get_read2()
		if len(primary_r2) == 0:
			return None
		elif len(primary_r2) == 1:
			return self._combine_short_CIGAR_elements(tol, primary_r2[0])	
		else:
			raise ValueError(f"found more than one primary alignment for {self[0].query_name} read 2")	

	def get_primary_not_r1_or_r2(self, tol):
		""" Return the primary read 1 alignment as a pysam.AlignedSegment object, or None if there isn't one """
		
		primary = self.get_primary().get_not_read1_or_read2()
		if len(primary) == 0:
			return None
		elif len(primary) == 1:
			return self._combine_short_CIGAR_elements(tol, primary[0])
		else:
			raise ValueError(f"found more than one primary alignment for {self[0].query_name} (not read 1 or read 2)")	
			
	def get_rearrangement_nm(self):
		""" 
		For an AlignmentPool consisting of different alignments of the same read, get an 
		edit distance for the case that those alignments cover the whole read
		
		First, we check each read to see if it has multiple mapped regions.  If it does,
		we split it into multiple reads, each of which has one of the mapped regions.
		
		Then remove any nested alignments using self.remove_query_nested().  Nested
		alignments are those that account for the same part of the read, or one accounts
		for a subset of the part of the read that is accounted for by the other.
		
		Finally, the rearrangement edit distnace is the sum of the edit distances for all 
		alignments, plus the number of bases in the read not accounted for by any alignment 
		(i.e. the number of bases in the gaps between those alignments) 
		
		Return a tuple of the edit distance and the processed (split, non-nested) AlignmentPool 
		(since it might be different to self after removing nested alignments)
		"""
		
		# check that we're dealing with the same query, and
		# check that all are read1, read2, or neither
		assert self._is_same_read()
		
		query_length = self[0].query_length
		
		alns = AlignmentPool()
		for aln in self:
			alns.append(self._copy_alignment(aln))
			
		alns._split_mapped()
		alns._remove_query_nested()
		
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
		nm += query_length - alns[-1].query_alignment_end
		
		return nm, alns
		
	def process_non_primary_alignments(self, tol, primary=None):
		"""
		add alignments from XA or SA of a primary alignment to the pool
		
		either provide a primary alignment, or do this based on the primary alignments
		already in the pool
		"""
		
		
		if primary is not None:
			assert self._is_same_read()
		
			self._process_non_primary_alignments(tol, primary)
			
			# remove any redundant alignments from pool
			# redundant alignments have same CIGAR, orientation (+/-) and pos		
			self.remove_redundant_alignments()
		
		else:
			for read in self:
				if read.is_secondary: 
					continue
				if read.is_supplementary:
					continue
				self._process_non_primary_alignments(tol, read)	
			
	def remove_redundant_alignments(self):
		""" 
		Remove redundant alignments (those that have same pos, CIGAR and strand) 
		Note that two alignments are also redundant if they have the same CIGAR, except that
		soft-clipped operations are replaced by hard-clipped operations
		"""
		
		assert self._is_same_read()
		
		if len(self) < 2:
			return
		
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
					i_found = [ind for ind in range(len(same)) if i in same[ind]]
					assert len(i_found) == 1
					ind = i_found[0]
					if j not in same[ind]:
						same[ind].append(j)
				elif j in redundant:
					j_found = [ind for ind in range(len(same)) if i in same[ind]]
					assert len(j_found) == 1
					if i not in same[ind]:
						same[ind].append(i)			
				else:
					same.append([i, j])
				
				redundant.add(i)
				redundant.add(j)
				
		# remove integrations from groups of redundant alignments, leaving only one
		# integration per group
		to_remove = []
		for group in same:
		
			# want to preferentially keep primary alignments, so if there's primary
			# alignments put them at the front of the list
			group = sorted(group, key=lambda x: self[x].is_supplementary or self[x].is_secondary)
			
			# on each loop, add one alignment to to_remove until there is only one left
			while len(group) > 1:
				removed = False
				
				for i in range(len(group)):
					
					# preferentially remove any alignments that have hard clips
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
			
	def _check_md_character(self, chr):
		"""
		Characters in MD tag can be integers, letters or '^'
		Return 'int' if inter, 'let' if letter, or 'del' if '^'
		"""
		if chr == '^':
			return 'del'
		
		try:
			int(chr)
			return 'int'
		except ValueError:
			return 'let'	
			
	def _combine_MD_with_CIGAR(self, cigartuples, md):
		"""
		Split a CIGAR and MD tag to get the CIGAR operations and their corresponding
		elements in the MD tag.  Combine into a list of tuples with length 3.
		In each tuple, first element is CIGAR op type, second element is CIGAR op length,
		and third element is list of MD elements
		For example CIGAR 1S249M and MD tag 10T2G7T5G221 should become
		[(4,1,[]), (0, 249, [10, 'T', 2, 'G', 7, 'T', 5, 'G', 221])]
		"""		
		
		if md is None:
			return [(op[0], op[1], []) for op in cigartuples]
		
		md_lst = self._split_md_tag(md)
		
		# if no md, just return empty lists for each element
		if len(md_lst) == 0:
			return [(op[0], op[1], []) for op in cigartuples]
		
		md_ind = 0
		md_cigartuples = []
		
		# loop over cigar operations and get corresponding 
		for op in cigartuples:
			
			cigar_op_MD = []
			# only ops that appear in MD tag are M and D
			
			# for a deletion, should be only one MD element for the CIGAR element
			if op[0] == 2:
				# deletion element looks like '^AGTG', where AGTG are deleted bases
				# number of deleted bases should match length of CIGAR op
				assert len(md_lst[md_ind]) - 1 == op[1]
				cigar_op_MD.append(md_lst[md_ind])
				md_ind += 1
			
			# for a mapped region, we need to collect the MD elements (matches (ints) and
			# single letter mismatches)
			elif op[0] == 0:
				md_len = 0
				while md_len < op[1]:
					md_op = md_lst[md_ind]
					
					try:
						# integer with a run of matched bases
						md_len += md_op
					except TypeError:
						# single letter mismatch
						assert len(md_op) == 1
						md_len += 1
										
					cigar_op_MD.append(md_op)
					md_ind += 1
				
				#it's possible to have md_len be longer than the operation, 
				# for example if two mapped regions are separated by an insertion
				if md_len > op[1]:
					cigar_op_MD[-1] -= md_len - op[1]
					md_ind -= 1
					md_lst[md_ind] -= cigar_op_MD[-1]

			
			md_cigartuples.append((op[0], op[1], cigar_op_MD))
		
		return md_cigartuples
				
	def _combine_short_CIGAR_elements(self, tol, read):
		"""
		Sometimes, short insertions, deletions and soft-clipped regions can cause a read
		to be missed because it doesn't meet our rather strict criteria for a chimeric read
		
		For example, viral CIGAR 1S100M49S and host CIGAR 101S49M would be missed as a
		potential integration because the viral cigar is clipped on both ends,
		even though it really looks like an integration.  
		
		To address this issue, combine any non-mapped short CIGAR elements ()
		that are shorter in length than a tolerance value (self.tol)
		
		Similar in spirit to '_simplify_CIGAR', but takes a slightly different approach (
		for example, deals differently with reads that are soft-clipped on both ends)
		
		We might later on want to split this alignment into smaller ones, for example in
		the case of a FullAlignment. When we do this, we will want an edit distance for
		each mapped region of the alignment. Therefore, when editing we should update the
		MD tag so that we can later compare query and reference for that mapped region.
		
 		We don't have all the information to do this properly (since the MD tag stores 
 		information about the reference, but we don't have the reference).  
 		We could solve this by asking the user for the reference as well and extracting 
 		the appropriate bases, but since we only care about the edit distance in this case, 
 		this is probably overkill.  
		
		So instead, just add 'X's to the MD tag when combining soft-clip elements with
		a mapped region.
		
		For deletion elements (D, 2), these are already included in MD tag but will be removed
		when we remove the element. Use the 'CO' (free-text comment tag) to keep track of
		how long these elements are and into which mapped region they were combined
		so that we can later add them into the edit distance if we split the read.
		For insertion elements (I, 1), these will be added to the MD tag as mismatches,
		but we also need to keep track of which mapped element they were combined into so 
		that we can later add them into the edit distance if we split the read, and to 
		account for their effect on coordinates in the reference.
		
		For insertions and deletions, use the CO tag to keep track of how many inserted
		and deleted bases were combined into a mapped region.  Do this in the format
		CO: '0:2:0,2:3:0,5:0:1' - a comma separated list of mapped regions in the new CIGAR, with
		the index of the mapped region in the cigartuples and the number of inserted and
		deleted bases merged separated by a colon.  For example, ''0:2:0,2:3:0,5:0:1' means
		that 2 inserted, 0 deleted bases were combined into the mapped region at index 0,
		3 inserted, 0 deleted were combined into the mapped region at index 2, 
		and 0 inserted, 1 deleted were combined into the mapped region at index 5.  Mapped 
		regions referred to in the CO should be AFTER
		editing the alignment - i.e. index 0 really should be a mapped region
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
				if curr_non_mapped > 0 and curr_non_mapped <= tol:
					to_delete.append((i, curr))
				curr_non_mapped = 0
				curr = []
				last_mapped = i
			# if not mapped
			else:
				curr.append(i)
				curr_non_mapped += read.cigartuples[i][1]
		
		# check last element
		if curr_non_mapped > 0 and curr_non_mapped <= tol:
			to_delete.append((last_mapped, curr))
		
		if len(to_delete) == 0:
			return read

		# next, combine short elements with nearest matched region
		
		# get MD tag for manipulating
		try:
			md = read.get_tag('MD')

			has_md = True
		except KeyError:
			has_md = False
			
		if has_md:
			tmp_cigartuples = self._combine_MD_with_CIGAR(read.cigartuples, md)
		else:
			tmp_cigartuples = [(i[0], i[1], []) for i in read.cigartuples]
				
		nm_offset = 0
		co = {}
		
		for elem in reversed(to_delete):

			i_matched = elem[0]
			assert read.cigartuples[i_matched][0] == 0
			assert tmp_cigartuples[i_matched][0] == 0
			i_del = elem[1]
			i_md = tmp_cigartuples[i_matched][2]
			
			# for co tag - keep track of number of I and D bases combined
			i_final_matched = i_matched
			try:
				num_ins, num_del = co[i_final_matched]
				del co[i_matched]
			except KeyError:
				num_ins, num_del = 0, 0
							
			# if cigar operation consumes query, need to add to matched region			
			matched = tmp_cigartuples[i_matched][1]
			for idx in reversed(i_del):
			
				## update cigar operation 
			
				# I (1) and S (4) operations consume query - add to matched length and edit distance offset
				if read.cigartuples[idx][0] in {1, 4}:
					matched += read.cigartuples[idx][1]
					
				# only soft-clipped operations need to be added to edit distance - 
				# insertions and deletions are already included
				if read.cigartuples[idx][0] == 4:
					nm_offset += read.cigartuples[idx][1]
					
				# if this is a deletion from reference (was previously in MD tag but will be 
				# removed), we will account for this in the CO tag so remove from NM offset
				if read.cigartuples[idx][0] == 2:
					nm_offset -= read.cigartuples[idx][1]
				
				## update MD tag
				# if we're combining a soft-clip or insertion with a matched region, 
				# add to MD tag
				if read.cigartuples[idx][0] in {4, 1}:
					
					# need to add bases to MD tag, but don't know what reference was at
					# this position.  So just add a string of 'X's, which need to be 
					# separated by 0's as per sam specification
					md_add = [0, 'X'] * read.cigartuples[idx][1] + [0]
					
					# if we're combining something after the mapped region
					if idx > i_matched:
					
						# don't need 0 at end
						md_add.pop(-1)
						# if last element of i_md is an int, don't need 0 at start
						if has_md:
							if isinstance(i_md[-1], int):
								md_add.pop(0)
					
						i_md = i_md + md_add
					
					# if we're combining something before the mapped region
					else:
						# don't need 0 at start
						md_add.pop(0)
						# if first element of i_md is an int, don't need 0 at end
						if has_md:
							if isinstance(i_md[0], int):
								md_add.pop(-1)
							
						i_md = md_add + i_md
					
				# for an insertion or deletion, keep track of the number of bases added
				# for adding to the CO tag
				if read.cigartuples[idx][0] in {1, 2}:
				
					if read.cigartuples[idx][0] == 1:
						num_ins += read.cigartuples[idx][1]
					else:
						num_del += read.cigartuples[idx][1]
					
				# if we're deleting something from before the mapped region, we also
				# need to offset our mapped region index in the CO tag
				if idx < i_matched:
					i_final_matched -= 1
				
			# update mapped length and md tag in cigartuples
			tmp_cigartuples[i_matched] = (0, matched, i_md)
			
			# remove elements
			for idx in reversed(i_del):
				tmp_cigartuples.pop(idx)

			# update CO tag
			co[i_final_matched] = (num_ins, num_del)
		
		# check that we don't now have two matched regions next to each other
		if len(tmp_cigartuples) > 1:
			for i in reversed(range(len(tmp_cigartuples)-1)):
				if tmp_cigartuples[i][0] == 0:
					if tmp_cigartuples[i+1][0] == 0:
						# add matched regions together
						matched = tmp_cigartuples[i][1] + tmp_cigartuples[i+1][1]
				
						# check that we wont't have two int md elements next to each other
						if has_md:
							if isinstance(tmp_cigartuples[i][2][-1], int):
								if isinstance(tmp_cigartuples[i+1][2][0], int):
									tmp_cigartuples[i+1][2][0] += tmp_cigartuples[i][2].pop(-1)
								
						# add MD for adjacent regions
						md =  tmp_cigartuples[i][2] + tmp_cigartuples[i+1][2]
						
						# we also need to adjust CO tag to make sure it refers to 
						# correct mapped region
						if i in co.keys() and i+1 in co.keys():
							# if we have both, combine them
							co[i] = (co[i+1][0] + co[i][0], co[i+1][1] + co[i][1])
							del co[i+1]
						
						# if we have only i+1 in co
						elif i+1 in co.keys():
							co[i] = co[i+1]
							del co[i+1]
						
						# remove one of regions
						tmp_cigartuples.pop(i+1)
						# update other with new matched and md
						tmp_cigartuples[i] = (0, matched, md)	
		
		# separate out md tag from cigartuples
		md = []
		for i in tmp_cigartuples:
			md += i[2]
		# check if we now have two numbers next to each other in the MD string 
		i = 1
		while i < len(md):
			if isinstance(md[i-1], int) and isinstance(md[i], int):
				md[i] = md[i] + md[i-1]
				md.pop(i-1)
				continue
			i += 1
		md = ''.join([str(i) for i in md])
		
		
		# if we didn't have an MD tag in the first place, don't add one
		if not has_md:
			md = None
		
		tmp_cigartuples = [(i[0], i[1]) for i in tmp_cigartuples]
		
		# join CO tag elements
		co = [f"{key}:{val[0]}:{val[1]}" for key, val in co.items()]
		co = ','.join(co)
			
		# if we're removing elements that consume the query and are before the first
		# mapped region for a forward read or after the last mapped region for a reverse 
		# read then we need to also adjust the mapping position
		deleted_inds = [tup[1] for tup in to_delete]
		#https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-a-list-of-lists
		deleted_inds = [item for sublist in deleted_inds for item in sublist]
		mapped = [ind for ind in range(len(read.cigartuples)) if read.cigartuples[ind][0] == 0]

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

		nm = read.get_tag('NM') + nm_offset
		
		return self._create_new_segment_from_primary(read, pos_offset=pos_offset, cigartuples=tmp_cigartuples,
														nm=nm, co=co, md=md)
			
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
		
	def _copy_alignment(self, read):
		"""
		Copy the attributes of a pysam.AlignedSegment into a new one
		"""

		return copy.deepcopy(read)	
		
	def _create_new_segment_from_primary(self, primary, ref_name=None, pos=None, 
											pos_offset=None, ori=None, 
											cigar=None, cigartuples=None, nm=None, 
											co=None, mapq=None, md=None):
		""" 
		Create a new pysam.AlignedSegment from an existing (primary) pysam.AlignnedSegment,
		but modify one or more of its properties
		"""
		
		# can't set both cigar (string) and cigartuples
		assert cigar is None or cigartuples is None
		
		a = self._copy_alignment(primary)
		
		if mapq is not None:
			a.mapping_quality = mapq
		
		if ori is None:
			ori = '-' if primary.is_reverse else '+'
		
		assert ori == '-' or ori == '+'
		
		# if primary alignment has different orientation to this one, 
		# need to reverse complement read sequence
		if primary.is_reverse != (ori == '-'):
			a.query_sequence = _reverse_complement(primary.query_sequence)
			a.is_reverse = not primary.is_reverse
#			a.query_qualities = copy.copy(primary.query_qualities)
#			a.query_qualities.reverse()


		if pos is not None:
			a.reference_start = pos

		assert a.reference_start >= 0
		if ref_name is None:
			assert a.reference_start <= primary.header.get_reference_length(primary.reference_name)
		else:
			assert a.reference_start <= primary.header.get_reference_length(ref_name)
		
		# use ZP tag to keep track of position offset as a result of manipulating read
		if pos_offset is not None:
			try:
				total_pos_offset = primary.get_tag('ZP') + pos_offset
			except KeyError:
				total_pos_offset = pos_offset
			if total_pos_offset != 0:
				a.set_tag('ZP', total_pos_offset)
		
		if ref_name is not None:
			a.reference_id = primary.header.get_tid(ref_name)
		
		if cigar is not None:
			a.cigartuples = self._convert_cigarstring_to_cigartuples(cigar)
		
		if cigartuples is not None:
			a.cigartuples = cigartuples
		
		if nm is not None:
			a.set_tag('NM', nm)
			
		if co is not None:
			a.set_tag('CO', co)
			
		if md is not None:
			a.set_tag('MD', md)

		# add OA tag
		orig_strand = '+' if not primary.is_reverse else '-'
		OA_update = (primary.reference_name, str(primary.reference_start), orig_strand, 
						primary.cigarstring, str(primary.mapping_quality), 
						str(primary.get_tag('NM')))
		OA_update = ",".join(OA_update) + ";"	
		try:
			a.set_tag('OA', primary.get_tag('OA') + OA_update)
		except KeyError:
			a.set_tag('OA', OA_update)

		return a	
		
	def _split_co(self, tag):
		"""
		Split co tag into a dict
		Tag has format : [{i_map}:{num_I}:{num_D}].+
		Split into format: {i_map: (num_i, num_D)}.+
		"""
		if tag is None:
			return None
		
		if tag == '':
			return None
			
		tag_dict = {}
		tag = tag.split(",")
		for elem in tag:
			elem = elem.split(":")
			assert len(elem) == 3
			tag_dict[int(elem[0])] = (int(elem[1]), int(elem[2]))
			
		return tag_dict
		
	def _get_edit_dist_for_mapped_region(self, read, map_op_idx):
		"""
		get the edit distance for a mapped region (at position map_op_idx in cigartuples)
		using the MD tag
		"""
		
		# check that this CIGAR operation exists and is mapped
		assert read.cigartuples[map_op_idx][0] == 0
		
		# if no MD tag, return None
		try:
			md = read.get_tag('MD')
		except KeyError:
			return None
			
		# combine cigartuples with MD tag to get part of MD tag for this element
		cigar_md = self._combine_MD_with_CIGAR(read.cigartuples, md)
		
		nm = 0
		for md_elem in cigar_md[map_op_idx][2]:
				
			# MD elements in a mapped region should only be single letter mismatches
			# or integers
			try:
				int(md_elem)
			except ValueError:
				assert len(md_elem) == 1
				nm += 1
		
		return nm		
		
	def _is_same_read(self):
		"""
		check if all the reads in self are the same read - this means
		that they all have the same query name, and read1/read2 status
		"""
		if len(self) == 1:
			return True
		
		if not all([self[0].query_name == read.query_name for read in self[1:]]):
			return False
			
		if not all([self[0].is_read1 == read.is_read1 for read in self[1:]]):
			return False

		if not all([self[0].is_read2 == read.is_read2 for read in self[1:]]):
			return False		
		
		return True

	def _process_non_primary_alignments(self, tol, primary):
		""" 
		Depending on aligner, secondary and supplementary alignments might be present in
		tags from the primary alignment, or as separate alignments, or both.  
		
		We might also want to simplify these alignments using self._combine_short_CIGAR_elements()
		"""
		
		self._process_SA_and_XA(primary)
		
		# combine any short elements in secondary alignments
		for i in range(len(self)):
			self[i] = self._combine_short_CIGAR_elements(tol, self[i])
	
	def _process_SA(self, SA, primary):
		""" 
		Create pysam.AlignedSegment objects for each XA alignment 
		
		SA has format (rname,pos,strand,CIGAR,mapQ,NM;)+
		"""
		
		SA_segs = AlignmentPool()
		
		SA = SA.split(";")
		SA = [i for i in SA if i != '']
		
		for aln in SA:
			self.add_read_from_SA(aln, primary)
				
	def _process_SA_and_XA(self, primary):
		""" Create AlignmentPool objects with alignments from XA and SA """
		
		# get supplementary alignments from SA tag
		try:
			SA = primary.get_tag('SA')
			self._process_SA(SA, primary)
		except KeyError:
			pass
		
		# get supplementary alignments from XA tag
		try:
			XA = primary.get_tag('XA')
			self._process_XA(XA, primary)
		except KeyError:
			pass				
	
	def _process_XA(self, XA, primary):
		""" 
		Create pysam.AlignedSegment objects for each XA alignment 
		
		XA is bwa-specific and has format (chr,pos,CIGAR,NM;)*
		"""
		
		XA = XA.split(";")
		XA = [i for i in XA if i != '']
		
		XA_segs = AlignmentPool()
		for aln in XA:
			self.add_read_from_XA(aln, primary)			
		
		
		return XA_segs
		
	def _remove_MD_tag(self, read):
		"""
		remove MD tag from read where it's not appropriate anymore
		"""
		# remove MD tag from read - don't apply to this alignment
		tags = read.get_tags()
		MD = [i for i in range(len(tags)) if tags[i][0] == 'MD']
		assert len(MD) == 0 or len(MD) == 1
		if len(MD) == 1:
			tags.pop(MD[0])
			
		read.tags = tags
		
		return read
		
	def _remove_query_nested(self):
		""" 
		Remove any alignments that are encompassed by other alignments.  That is,
		remove any alignments that cover the same part of the query, keeping the alignment 
		that covers the larger part of the query
		
		For example, considering two alignments, one 150M and 130M20S, retain only 150M.
		
		For two alignments, 50M100S and 40S110M, keep both because neither alignment
		completely encompasses the other
		
		If two alignments are equivalent (i.e. cover same part of read), remove the one
		with the higher edit distance
		
		Note that this changes the order of alignments in the pool (they will be sorted
		by query alignment start position)
		
		It also reverses any reads that are mapped in the reverse orientation, which will 
		mess up the mapping positions
		"""
		assert self._is_same_read()
		
		if len(self) < 2:
			return
			
		assert all([read.query_name == self[0].query_name for read in self[1:]])
		
		# reverse CIGAR for any reads that are mapped in the reverse orientation
		for i in range(len(self)):
			
			if self[i].is_reverse:
				self[i] = self._reverse_cigar(self[i])

		
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
						if self[i].get_tag('NM') > self[i-1].get_tag('NM'):
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
		
	def _reverse_cigar(self, read):
		"""
		Create a new pysam.AlignedSegment with the CIGAR string reversed.  Nothing
		else is changed, so this alignment doesn't really make sense anymore. Use only
		for finding edit distance for rearrangement!
		"""	
		cig = read.cigartuples
		cig.reverse()
		
		a = pysam.AlignedSegment()
		a.cigartuples = cig
		a.set_tag('NM', read.get_tag('NM'))
		
		return a
		
	def _reverse_co(self, tag, cigar_len):
		"""
		reverse a CO tag, which has format
		"{i_map}:{num_I}:{num_D}"
		"""
		# we need to know the length of cigartuples 
		
		tag = tag.split(",")
		if tag != ['']:
			for i in range(len(tag)):
				
				map_i, I, D = tag[i].split(":")
				
				# get index in reversed cigar
				map_i = abs(cigar_len - int(map_i))
				
				tag[i] = f"{map_i}:{I}:{D}"
				
			tag = tag.join(",")
			
		return tag	
		
	def _split_mapped(self):
		"""
		For each alignment in the pool, check if it has multiple mapped regions. 
		If it does, remove the alignment from the pool, create multiple reads which each
		have one of the mapped regions, and add these to the pool.
		"""
		
		for i, read in enumerate(self):
			mapped = [op for op in read.cigartuples if op[0] == 0]
			if len(mapped) <= 1:
				continue
			
			self.pop(i)
			
			split_reads = self._split_read(read)

			for splt in split_reads:
				self.insert(i, splt)		
				
	def _split_md_tag(self, md):
		"""
		Split MD tag into individual elements for ease of manipulation
		e.g. '10A5^AC6' should become 
		[10, 'A', 5, '^AC', 6]
		"""

		md_lst = []
		
		# if no md tag, return empty list
		if len(md) == 0:
			return md_lst
			
		val = md[0]
		# keep track of what the previous character was
		prev_type = self._check_md_character(md[0])

		for i in md[1:]:
			curr_type = self._check_md_character(i)
			# if previous type was a number
			if prev_type == 'int':
				# if this is also a number, it's part of the same number
				if curr_type == 'int':
					val = val + i
				# if it's not part of a number add the previous number and start the new 
				# element
				else:
					md_lst.append(int(val))
					val = i
			# if this is part of a deletion
			elif prev_type == 'del':
				# if it's a letter, then add to the deletion
				if curr_type == 'let':
					val = val + i
					curr_type = 'del'
				# otherwise has to be a number - add the deletion to the list and start
				# the new number
				else:
					md_lst.append(val)
					val = i
			# otherwise has to be a single letter mismatch
			else:
				md_lst.append(val)
				val = i
			
			prev_type = curr_type
		
		# append last element
		try:
			md_lst.append(int(val))
		except ValueError:
			md_lst.append(val)
				
		return md_lst	
		
	def _split_read(self, read):
		"""
		Split a read with multiple mapped regions into multiple alignments, one for
		each mapped region
		
		Each mapped region will be soft-clipped either side, regardless of what was
		previously in the CIGAR
		
		e.g. 20S50M30I30M20D30M20S will become three alignments:
		20S50M110S
		100S30M50M
		130S30M20S
		"""	
		
		mapped_idxs = [i for i in range(len(read.cigartuples)) if read.cigartuples[i][0] == 0]
		
		assert len(mapped_idxs) > 1
		
		new_reads = AlignmentPool()
		
		# get CO tags, if it exists
		try:
			co = read.get_tag('CO')
			co = self._split_co(co)
		except KeyError:
			co = {}
			
		try:
			md = read.get_tag('MD')
		except KeyEror:
			md = None

		pos_offset = 0
		for i, cigar_op_idx in enumerate(mapped_idxs):
			
			# create new cigartuples
			new_cigartuples = []
			# get the number of bases consumed by CIGAR operations before this mapped region
			num_before = [i[1] for i in read.cigartuples[:cigar_op_idx] if i[0] in (0, 1, 4)]
			num_before = sum(num_before)
			
			# soft-clipped region preceeds mapped region
			if num_before != 0:
				new_cigartuples.append((4, num_before))
			
			# add mapped region
			mapped = read.cigartuples[cigar_op_idx][1]
			new_cigartuples.append((0, mapped))
		
			# get number of bases after
			num_after = len(read.query_sequence) - mapped - num_before
			
			if num_after > 0:
				new_cigartuples.append((4, num_after))
				
			# get an edit distance for this M operation only
			map_nm = self._get_edit_dist_for_mapped_region(read, cigar_op_idx)

			# get mapping position for this M operation
			try:
				pos = read.get_blocks()[i][0] + pos_offset
			except OverflowError:
				pdb.set_trace()
				
			# If we combined inserted or deleted elements into this mapped region, we 
			# need to offset the mapping position for all elements following that mapped
			# region
			try:
				pos_offset -= co[cigar_op_idx][0]
				pos_offset += co[cigar_op_idx][1]
			except KeyError:
				pass
				
			# we also want to replace the CO tags for the original read
			# with entries appropriate for just this mapped element
			try:
				op_co = f"{0 if num_before == 0 else 1}:{co[cigar_op_idx][0]}:{co[cigar_op_idx][1]}"
			except KeyError:
				if co is None or co == {}:
					op_co = None
				else:
					op_co = ''
				
			# get MD tag for just this mapped region, if there was one
			combined_cigar_md = tmp_cigartuples = self._combine_MD_with_CIGAR(read.cigartuples, md)
			op_md = combined_cigar_md[cigar_op_idx][2]
			if len(op_md) > 0:
				op_md = "".join([str(i) for i in op_md])
			else:
				op_md = None
			
			a = self._create_new_segment_from_primary(read, pos=pos, 
														cigartuples=new_cigartuples, 
														nm=map_nm, co=op_co, md=op_md)
			
			new_reads.append(a)

		return new_reads	

					
class ChimericIntegration:
	""" 
	A class to store/calculate the properties of a simple chimeric integration 
	(ie mapped/clipped and clipped/mapped) 
	"""
	__slots__ = 'map_thresh', 'primary', 'hread', 'vread', 'hsec', 'vsec', 'chr', 'virus', 'verbose', 'tol', 'nm_pc', 'nm_diff'
	
	def __init__(self, host, virus, host_sec = AlignmentPool(), virus_sec = AlignmentPool(),  
					 tol = 3, map_thresh=20, nm_diff=None,
					  nm_pc=None, primary=True, verbose=False):
		""" 
		Given a host and virus read that are chimeric, calculate the properties of the integration.
		
		If the ChimericIntegration is part of a FullIntegration, 'side' should be 'left' 
		or 'right', otherwise None if the integration encompasses the whole read.
		"""

		self.map_thresh = map_thresh
		self.tol = tol
		self.nm_pc = nm_pc
		self.nm_diff = nm_diff
		
		assert self.map_thresh > 0
		assert self.tol >= 0
		
		if self.nm_pc is not None:
			assert self.nm_pc >= 0 and self.nm_pc <= 1
		if self.nm_diff is not None:
			assert self.nm_diff >= 0
		
		self.primary = primary
		self.verbose = verbose

		self.hread = host

		self.vread = virus

		# simple chimeric read has two CIGAR elements
		assert self._is_chimeric()
		
		self.hsec = host_sec	
		self.vsec = virus_sec
		
		self.chr = host.reference_name
		self.virus = virus.reference_name
					
	def get_ambig_seq(self):
		""" Get sequence of ambiguous base(s) from read """

		ambig = self._get_ambig_coords_read()
		
		if not self.hread.is_reverse:
			return self.hread.query_sequence[ambig[0]:ambig[1]]
		else:
			return _reverse_complement(self.hread.query_sequence[ambig[0]:ambig[1]])

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
		
		coords = self._get_ambig_coords_ref(self.hread)
		
		return coords

	def get_host_edit_dist(self):
		""" Get the edit distance for the host alignment """
		
		return self._get_edit_dist(self.hread)

	def get_host_mapq(self):
		""" Get mapping quality for host alignment """
		return self.hread.mapping_quality

	def get_host_seq(self):
		""" 
		Get the part of the read aligned to the host 
		(including ambiguous bases if overlap, excluding if gap)
		"""
		
		if self.hread.is_reverse:
			return _reverse_complement(self.hread.query_alignment_sequence)
		else:
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
		AltLocs: A list of the alternative integrations
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
			'AltLocs': ";".join([integ.short_str() for integ in alt_ints]),
			'Type': self.get_integration_type(),
			'HostMapQ': self.get_host_mapq(),
			'ViralMapQ': self.get_viral_mapq(),
			'ReadID': self.get_read_id(),
			'ReadSeq': self.get_read_sequence()
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
			
	def get_read_sequence(self):
		"""
		Get the sequence of the read
		"""
		if not self.hread.is_reverse:
			return self.hread.query_sequence
		else:
			return _reverse_complement(self.hread.query_sequence)
	
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
		
		coords = self._get_ambig_coords_ref(self.vread)	
		
		return coords
		
	def get_viral_edit_dist(self):
		""" Get the edit distance for the viral alignment """
		
		return self._get_edit_dist(self.vread)
	
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
		""" 
		Get the part of the read aligned to the host 
		(including ambiguous bases if overlap, excluding if gap)
		"""
		
		if self.vread.is_reverse:
			return _reverse_complement(self.vread.query_alignment_sequence)
		else:
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
		
		Normally host and virus reads should be the same length, because they are alignments
		of the same read.  However, if this ChimericIntegration comes from a FullIntegration,
		we have split up each host and virus read into two parts, which may be of different
		lengths.  In this case, we need to allow for reads of different lengths 
		
		"""

		if self.get_integration_orientation() == 'hv':

			hop = self.hread.cigartuples[0]
			if self.vread.is_reverse == self.hread.is_reverse:
				vop = self.vread.cigartuples[0]
			else:
				vop = self.vread.cigartuples[-1]
				
		else:
			hop = self.hread.cigartuples[-1]
			if self.vread.is_reverse == self.hread.is_reverse:
				vop = self.vread.cigartuples[-1]
			else:
				vop = self.vread.cigartuples[0]				
		
		# subtract mapped from soft-clipped 
		if hop[0] == 0:
			assert vop[0] == 4
			ambig = hop[1] - vop[1]
		else:
			assert hop[0] == 4
			ambig = vop[1] - hop[1]
			
		if absolute:
			return abs(ambig)
		
		return ambig
		
	def is_host_ambig_loc(self):
		""" 
		Check if an integration has ambiguous location in host genome (i.e. host part of
		read has multiple equivalent alignments to different parts of host genome)
		"""
		return self._is_ambiguous_location(self.hread, self.hsec)
	
	def is_possible_translocation(self):
		""" 
		A read is a possible host translocation if more bases in the read can be 
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
		
	def set_sec_alns(self, hsec, vsec):
		"""
		if initialised without secondary or supplementary alignments, add them
		"""
		
		hsec.process_non_primary_alignments(self.tol, self.hread)
		vsec.process_non_primary_alignments(self.tol, self.vread)
		
		self.hsec = hsec
		self.vsec = vsec
		
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
		
	def _filter_alt_ints(self, alt_ints):
		"""
		Filter alternate integrations using self.alt_diff and self.alt_pc
		"""

		# only keep alternative integration sites that are close in edit distance
		# to primary
		if self.nm_pc is not None:
			alt_ints = [int for int in alt_ints if int.get_total_edit_dist() == 0 or self.get_total_edit_dist()/int.get_total_edit_dist() >= self.nm_pc]
		
		if self.nm_diff is not None:
			alt_ints = [int for int in alt_ints if int.get_total_edit_dist() - self.get_total_edit_dist() <= self.nm_diff]
					
		return alt_ints
		
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
		
		if len(self.hsec) == 0 and len(self.vsec) == 0:
			return alt_ints
		
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
					alt_int = ChimericIntegration(hst, vrs, tol = self.tol,
									map_thresh = self.map_thresh, primary=False)
					alt_ints.append(alt_int)
					
				except AssertionError:
					pass

		alt_ints = self._swap_self_with_primary(alt_ints, host_alns, virus_alns)
		
		alt_ints = self._filter_alt_ints(alt_ints)
	
		return alt_ints

	def _get_ambig_coords_read(self):
		""" Get coordinates of ambiguous sequence relative to read """
		
		# if host/virus, want end of host alignment if forward
		if self.get_integration_orientation() == 'hv':
			return sorted([
				self.hread.query_alignment_end,
				self.hread.query_alignment_end - self.num_ambig_bases()
				])
		# if virus/host, want start of host alignment if forward
		else:
			return sorted([
				self.hread.query_alignment_start,
				self.hread.query_alignment_start + self.num_ambig_bases()
				])	
			
	def _get_ambig_coords_ref(self, read):
		""" 
		Return a tuple of (start, stop) of ambiguous bases in host or virus 
		(i.e. in coords in reference) 
		"""
		# decide if aligned part of read is first or second
		if read.cigartuples[0][0] == 0:
			first = True
		else:
			first = False

		# if read has ZP tag, use the position offset
		try:
			pos_offset = read.get_tag('ZP')
		except KeyError:
			pos_offset = 0
			
		# if read has deleted or inserted bases, we need to offset the end position
		try:
			co = read.get_tag('CO')
		except KeyError:
			co = ''
		end_pos_offset = self._get_co_D_bases(co) - self._get_co_I_bases(co)
			
		# get start and end of aligned part of read in host coords
		aligned_coords = (read.get_blocks()[0][0] + pos_offset, 
							read.get_blocks()[-1][-1] + pos_offset + end_pos_offset)

		# if aligned part of read is first
		if first:
			start = aligned_coords[1]
			stop = aligned_coords[1]
			
			if self.gap_or_overlap() == 'gap':
				stop += self.num_ambig_bases(absolute=True)
			elif self.gap_or_overlap() == 'overlap':
				start -= self.num_ambig_bases(absolute=True)

		# if aligned part of read is second
		else:
			start = aligned_coords[0]
			stop = aligned_coords[0]
			
			if self.gap_or_overlap() == 'gap':
				start -= self.num_ambig_bases(absolute=True) 
			elif self.gap_or_overlap() == 'overlap':
				stop += self.num_ambig_bases(absolute=True)

		# check that the start and stop are within the bounds of the chromosome/reference
		if self.gap_or_overlap() != 'gap':
			# for a clean junction or overlap, we should definitely be in the reference
			try:
				assert start >= 0
				assert stop <= read.header.get_reference_length(read.reference_name)
			except AssertionError:
				pdb.set_trace()
		else:
			# if a gap, we may have gone beyond reference - adjust coordinates
			if start < 0:
				start = 0
			if stop > read.header.get_reference_length(read.reference_name):
				stop = read.header.get_reference_length(read.reference_name)
				
		return start, stop
		
	def _get_co_bases(self, read):
		"""
		get total number of bases in CO tag
		"""	
		try:
			co = read.get_tag('CO')
		except KeyError:
			return 0
			
		if co == '':
			return 0

		return self._get_co_D_bases(co) + self._get_co_I_bases(co)

	def _get_co_D_bases(self, tag):
		"""
		get total number of deleted bases in CO tag
		"""
		dele = 0
		if tag == '':
			return dele
		tag = tag.split(",")
		for op in tag:
			dele += int(op.split(":")[2])
			
		return dele
		
	def _get_co_I_bases(self, tag):
		"""
		get total number of inserted bases in CO tag
		"""
		ins = 0
		if tag == '':
			return ins
		tag = tag.split(",")
		for op in tag:
			ins += int(op.split(":")[1])
			
		return ins
		
	def _get_edit_dist(self, read):
		"""
		Get edit distance for a read
		"""
			
		try:
			nm = read.get_tag("NM")
		except KeyError:
			return None
			
		# deleted bases combined into the mapped region are not accounted for in NM
		try:
			co = read.get_tag("CO")
		except KeyError:
			return nm
		
		return nm + self._get_co_D_bases(co)
				
	def _is_ambiguous_location(self, read, sec):
		"""  
		An alignment is ambiguous if there's a secondary or supplementary alignment
		that is equivalent to the primary one.
		An equivalent alignment is one with the same CIGAR and edit distance
		"""	

		if len(sec) == 0:
			return False
			
		# get edit distance for primary alignment
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
		
		# if there's only one alignment in pool after removing nested alignments,
		# can't be rearrangement
		if len(all) == 1:
			return False
			
		# compare integration edit distance with rearrangement edit distance
		int_nm = self.get_total_edit_dist()
		
		# if we have a nm_diff or nm_pc, this is our threshold for how much more int_nm
		# has to be compared to rearrange_nm to say that this is not a rearrangement
		if self.nm_diff is not None:
			return rearrange_nm - int_nm <= self.nm_diff
		elif self.nm_pc is not None:
			if rearrange_nm != 0:
				return int_nm / rearrange_nm >= self.nm_pc
			# we get a divide by zero error if rearrange_nm is zero, but if this is the case
			# then it's probably a rearrangement
			else:
				return True
		else:
			return rearrange_nm >= int_nm
					
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
			
		# soft-clipped regions must be at least 20bp
		if (self.hread.query_length  - self.hread.query_alignment_length) < self.map_thresh:
			return False
		if (self.vread.query_length  - self.vread.query_alignment_length) < self.map_thresh:
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
						
	def _swap_self_with_primary(self, alt_ints, host_alns, virus_alns):
		"""
		Check each alternate integration to see if it has an edit distance lower than 
		the primary edit distance.  If it does, the one with the lower edit distance
		becomes the primary and the other becomes an alternate
		"""
		# if this integration is a primary one, but we found an alternative integration with
		# a lower edit distance, then the alternate one becomes primary
		
		
		primary_nm = self.get_total_edit_dist()
		if self.primary:
		
			for i in range(len(alt_ints)):
				
				alt_nm = alt_ints[i].get_total_edit_dist()
			
				if alt_nm < primary_nm:
			
					if self.verbose:
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
	__slots__ = 'hread1', 'hread2', 'vread1', 'vread2', 'hsec1', 'hsec2', 'vsec1', 'vsec2', 'tlen'
	def __init__(self, host_r1, host_r2, virus_r1, virus_r2,
					host_sec1 = AlignmentPool(), host_sec2=AlignmentPool(), 
					virus_sec1=AlignmentPool(), virus_sec2=AlignmentPool(),
					map_thresh=20, tlen=0, tol=3, nm_diff=None,
					  nm_pc=None, primary=True, verbose=False):
					
		assert tlen >= 0
		assert tol >= 0
		assert map_thresh >= 0
		if nm_pc is not None:
			assert nm_pc >=0 and nm_pc <= 1
		if nm_diff is not None:
			assert nm_diff >= 0
		
		self.map_thresh = map_thresh
		self.tlen = int(tlen)
		self.tol = int(tol)
		self.nm_pc = nm_pc
		self.nm_diff = nm_diff
		
		self.primary = primary
		self.verbose = verbose
		
		self.hread1 = host_r1
		self.hread2 = host_r2
		self.vread1 = virus_r1
		self.vread2 = virus_r2
		
		self.hsec1 = host_sec1
		self.hsec2 = host_sec2
		self.vsec1 = virus_sec1
		self.vsec2 = virus_sec2	

		
		# note that self._is_discordant() also assigns self.hread and self.vread
		assert self._is_discordant(host_sec1, host_sec2, virus_sec1, virus_sec2,
									host_r1, host_r2, virus_r1, virus_r2)
		
		self.get_host_chr()
		self.get_viral_ref()	
		
	def gap_or_overlap(self):
		""" For a discordant pair, OverlapType is always discordant """
		return 'discordant'
	
	def get_ambig_seq(self):
		""" No ambiguous bases for a discordant pair"""

		return None

	def get_coords(self, read):
		"""
		Get the likely coordinates of integration.  This is the end of the aliged read
		that faces the other read, plus the mean insert size
		"""

		# mean insert length is mean template length minus sum of R1 and R2 lengths
		insert = self.tlen - self.vread.query_length - self.hread.query_length
		# but can't be less than zero
		if insert < 0:
			insert = 0
			
		# if we modified the read by combining small CIGAR elements, we modified
		# the mapping positions, in which case we need to apply an offet
		try:
			pos_offset = read.get_tag('ZP')
		except KeyError:
			pos_offset = 0
			
		# if read has deleted or inserted bases, we need to offset the end position
		try:
			co = read.get_tag('CO')
		except KeyError:
			co = ''
		end_pos_offset = self._get_co_D_bases(co) - self._get_co_I_bases(co)
			
		# get start and end of aligned part of read in host coords
		try:
			blocks = (read.get_blocks()[0][0] + pos_offset, 
							read.get_blocks()[-1][-1] + pos_offset + end_pos_offset)
		except:
			print(read.query_name)
		# we want the coordinate of the host/virus base closest to integration
		if read.is_reverse:
			stop = blocks[0]
			start = stop - insert
		else:
			start = blocks[-1]
			stop = start + insert
			
		# check the start isn't less than zero
		if start < 0:
			start = 0
		if stop < 0:
			stop = 0

		# check that the stop isn't more than the length of the reference
		reference_length = read.header.get_reference_length(read.reference_name)
		if stop > reference_length:
			stop = reference_length
			
		if start > reference_length:
			start = reference_length
		
		return (start, stop)

	def get_host_coords(self):
		""" Get likely coordinates of integration in host chromosome """
		
		return self.get_coords(self.hread)

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
		
	def get_read_sequence(self):
		"""
		Get the sequence of the reads - combine R1 and R2 with xxx in the middle
		"""
		return self.hread.query_sequence  + ";" + self.vread.query_sequence
		
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
    	or this:
    	
		              read 1                 read 2
		host:  xxxxxxxxxxxxxxxxxxxxx xxxoooooooooooooooooo
		virus: ooooooooooooooooooooo oooxxxxxxxxxxxxxxxxxx
                                     ^^^    	
    	
        In this case, these soft-clipped bases (top) and mapped bases (bottom) 
        are probably not unexpected - they indicate that the host/virus junction lies just 
        before the end of read 1.  In practice, we are not likely to actually see evidence 
        of arrangement though, because there is a minimum length of alignments (e.g. the 
        seed length of the aligner).  For example, if the junction is 5bp from the right of 
        read 1, all we will see is that the viral alignment is soft-clipped (5bp), but since 
        5bp is smaller than our seed length we won't see the 5bp alignment in the host
        alignment for read 1.  So we might conclude that we should add 5 to the edit 
        distance for the 5 soft-clipped bases on the right end of the viral alignment for 
        read 1, but we'd probably be wrong.
        
        We therefore don't have enough information to calculate an edit distance analogous
        to the chimeric case for a discordant integration (most of the time). We could
        output the sum of the host and viral alignment edit distances, but I think this
        is misleading because it doesn't (and can't) capture the whole picture. But,
        a total edit distance would be useful for deciding which of a list of alternate
        integrations (constructed using secondary alignments) are more likely.  So 
        output a
		"""
		
		# check for unmapped bases on the part of the read that is furtherest from the
		# integration
		h1_map = self._is_mapped(self.hread1)
		
		# first, add the edit distances of the two mapped reads
		# this will include insertions combined 
		nm = self.get_viral_edit_dist() + self.get_host_edit_dist()
		
		# check if there are any bases mapped at one end of the unmapped read (closest
		# to the mapped read) that are mapped in the other alignment
		if h1_map:
		
			nm += self._check_overlap_in_unmapped(self.hread1, self.vread1)
			nm += self._check_overlap_in_unmapped(self.vread2, self.hread2)
		else:
			nm += self._check_overlap_in_unmapped(self.vread1, self.hread1)
			nm += self._check_overlap_in_unmapped(self.hread2, self.vread2)
		
		return nm
	
	def get_viral_coords(self):
		""" Get likely coordinates of integration in viral reference """		

		return self.get_coords(self.vread)
		
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
			
	def set_sec_alns(self, hsec1, hsec2, vsec1, vsec2):
		"""
		if initialised without secondary or supplementary alignments, add them
		"""
		
		hsec1.process_non_primary_alignments(self.tol, self.hread1)
		hsec2.process_non_primary_alignments(self.tol, self.hread2)
		vsec1.process_non_primary_alignments(self.tol, self.vread1)
		vsec2.process_non_primary_alignments(self.tol, self.vread2)
		
		self.hsec1 = hsec1
		self.hsec2 = hsec2
		self.vsec1 = vsec1
		self.vsec2 = vsec2
		
		if self._is_mapped(self.hread1):
			self.hsec = self.hsec1
			self.vsec = self.vsec2
		else:
			self.hsec = self.hsec2
			self.vsec = self.vsec1
			
	def _check_overlap_in_unmapped(self, mapped, unmapped):
		"""
		Check if there are any soft-clipped bases in the 'mapped' read that actually overlap
		with mapped bases in the 'unmapped' one, e.g.
		
		              read 1                 read 2
		host:  xxxxxxxxxxxxxxxxxxxxx xxxoooooooooooooooooo
		virus: ooooooooooooooooooooo oooxxxxxxxxxxxxxxxxxx		
									 ^^^
									 
		Where o is mapped and x is unmapped.
		
		The following is also fine - indicates overlap in read 2
		
		              read 1                 read 2
		host:  xxxxxxxxxxxxxxxxxxxxx ooooooooooooooooooooo
		virus: ooooooooooooooooooooo oooxxxxxxxxxxxxxxxxxx	
									 ^^^

		The following is not fine - soft-clipped bases in R2 in host are not mapped in virus
		
		              read 1                 read 2
		host:  xxxxxxxxxxxxxxxxxxxxx xxxoooooooooooooooooo
		virus: ooooooooooooooooooooo xxxxxxxxxxxxxxxxxxxxx	
									 ^^^		
		Return the number of bases in this read that should be added to the edit distance.
		This is any mapped bases in the unmapped read that don't overlap in the read
		with soft-clipped bases in the mapped read, as well as any soft-clipped bases in
		the mapped read that don't overlap with soft-clipped bases in the overlapped read.
		In other words, it's any bases that don't match the first two scenarios above.
		"""
		
		# get index soft-clipped elements in mapped read
		soft_map = [ind for ind in range(len(mapped.cigartuples)) if mapped.cigartuples[ind][0] == 4]
			
		# get index of mapped elements in unmapped read
		try:
			map_unmap = [ind for ind in range(len(unmapped.cigartuples)) if unmapped.cigartuples[ind][0] == 0]
		# if read is completely unmapped, cigartuples is None
		except TypeError:
			map_unmap = []
		
		# if there aren't any of either, return 0
		if len(map_unmap) == 0 and len(soft_map) == 0:
			return 0

		# if we don't have any mapped bases in the unmapped read, then any soft-clipped
		# bases in the mapped read are unexplained
		if len(map_unmap) == 0:
			total_mapped = mapped.query_alignment_end - mapped.query_alignment_start
			return mapped.infer_query_length() - total_mapped
				
		# check if the mapped part of the unmapped read is in the 'correct' part - i.e.
		# closest part of read to integration junction / other read
		if unmapped.is_reverse:
			# if mapped in reverse, mapped should be first element in cigar
			if map_unmap[0] == 0:

				mapped_correct_pos = True
				
				# we won't consider these bases when calculating edit distance
				map_unmap.pop(0)
				
				# if there are soft-clipped bases in the same part of the mapped read,
				# we won't consider these either
				if len(soft_map) > 0:
					if mapped.is_reverse:
						if soft_map[-1] == len(mapped.cigartuples) - 1:
							soft_map.pop(-1)
					else:
						if soft_map[0] == 0:
							soft_map.pop(0)			
				
			else:
				mapped_correct_pos = False
		else:
			# otherwise, mapped should be last element in cigar
			if map_unmap[-1] == len(unmapped.cigartuples) - 1:
				
				mapped_correct_pos = True	
				
				map_unmap.pop(-1)	
				
				# also check if there's a matching soft-clipped part of read
				if len(soft_map) > 0:
					if mapped.is_reverse:
						if soft_map[0] == 0:
							soft_map.pop(0)
					else:
						if soft_map[-1] == len(mapped.cigartuples) - 1:
							soft_map.pop(-1)	
			else:
				mapped_correct_pos = False		
		
		# if mapped part of unmapped read isn't
		# in correct place, then we just add all unmapped bases in mapped read and
		# mapped bases in unmapped read
		if not mapped_correct_pos:
			# get number of mapped bases in unmapped read
			nm = unmapped.query_alignment_end - unmapped.query_alignment_start
			
			# get number of unmapped bases in mapped read
			total_mapped = mapped.query_alignment_end - mapped.query_alignment_start
			nm += mapped.infer_query_length() - total_mapped
			
			return nm
			
		# otherwise, count bases in remaining soft-clips and mapped regions
		nm = 0
		for idx in soft_map:
			nm += mapped.cigartuples[idx][1]
		for idx in map_unmap:
			nm += unmapped.cigartuples[idx][1]
			
		return nm
					
	def _is_discordant(self, host_sec1, host_sec2, virus_sec1, virus_sec2, 
							orig_hread1, orig_hread2, orig_vread1, orig_vread2):
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

		self.hsec1 = host_sec1
		self.hsec2 = host_sec2
		self.vsec1 = virus_sec1
		self.vsec2 = virus_sec2	
			
		# assign self.hread and self.vread to help with some methods inherited from ChimericIntegration
		if h1_map:
			
			self.hread = self.hread1
			self.vread = self.vread2
			self.hsec = self.hsec1
			self.vsec = self.vsec2			
		else:
			self.hread = self.hread2
			self.vread = self.vread1
			self.hsec = self.hsec2
			self.vsec = self.vsec1		

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
		
	def _is_possible_rearrangement(self, primary, sec):
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
		
		haln = AlignmentPool(self.hsec)
		valn = AlignmentPool(self.vsec)
		
		if hread1_map:
			
			haln.append(self.hread1)
			haln.remove_redundant_alignments()
			
			valn.append(self.vread2)
			valn.remove_redundant_alignments()
			
			h2 = self.hread2
			v1 = self.vread1
		
		else:

			haln.append(self.hread2)
			haln.remove_redundant_alignments()

			valn.append(self.vread1)
			valn.remove_redundant_alignments()
			
			h1 = self.hread1
			v2 = self.vread2

		if len(haln) == 1 and len(valn) == 1:
			return alt_ints
		
		for h_i, v_i in itertools.product(range(len(haln)), range(len(valn))):
				
				# don't do primary vs primary
				if haln[h_i] == self.hread and valn[v_i] == self.vread:
					continue
					
				if hread1_map:
					h1 = haln[h_i]
					v2 = valn[v_i]
				else:
					h2 = haln[h_i]
					v1 = valn[v_i]
				
				# if host and virus reads constitute a valid integration, add to list
				try:
					alt_int = DiscordantIntegration(h1, h2, v1, v2, 
										map_thresh = self.map_thresh, primary = False)
	
				except AssertionError:
					alt_int = None
			
				if alt_int is not None:
					# replace host and virus alignments in haln/valn with the edited version
					haln.pop(h_i)
					valn.pop(v_i)
					if hread1_map:
						haln.insert(h_i, h1)
						valn.insert(v_i, v2)
					else:
						haln.insert(h_i, h2)
						valn.insert(v_i, v1)						
				
					alt_ints.append(alt_int)				

		alt_ints = self._swap_self_with_primary(alt_ints, haln, valn)
		
		alt_ints = self._filter_alt_ints(alt_ints)
			
		return alt_ints
		
	def _get_last_mapped(self, read):
		"""
		Get the number of bases in the query between the last mapped region and the end
		of the read.  Consider only cigar operations that consume query (I/S)
		"""
		cigar = read.cigartuples
		mapped = [ind for ind in range(len(cigar)) if cigar[ind][0] == 0]
		last_mapped = mapped[-1]
		
		return sum([op[1] for op in cigar[last_mapped:] if op[0] in (1, 4)])
		
	def _get_first_mapped(self, read):
		"""
		Get the number of bases in the query between the start of the read and the first
		mapped region.  Consider only cigar operations that consume query (I/S)
		"""		
		cigar = read.cigartuples
		mapped = [ind for ind in range(len(cigar)) if cigar[ind][0] == 0]
		first_mapped = mapped[0]
		
		return sum([op[1] for op in cigar[first_mapped:] if op[0] in (1, 4)])		

class FullIntegration(ChimericIntegration):
	"""
	A 'full integration' is a read that looks like it contains both junctions (left and 
	right) of an integration.
	"""
	__slots__ = 'left', 'right'
	def __init__(self, host, virus, host_sec = AlignmentPool(), virus_sec = AlignmentPool(), 
					map_thresh=20, tol=3, nm_diff=None,
					  nm_pc=None, primary=True, verbose=False):
				
		assert map_thresh > 0
		assert tol >= 0
		if nm_pc is not None:
			assert nm_pc >=0 and nm_pc <= 1
		if nm_diff is not None:
			assert nm_diff >= 0
					
		self.map_thresh = map_thresh
		self.tol = tol	
		self.nm_pc = nm_pc
		self.nm_diff = nm_diff
		
		self.primary = primary	
		
		self.verbose = verbose
		
		self.hread = host
		self.vread = virus
		
		assert self._is_chimeric()
		
		self.hsec = host_sec
		self.vsec = virus_sec
		
		self.chr = host.reference_name
		self.virus = virus.reference_name	
		
		self._assign_junctions()
		
	def _assign_junctions(self):
		"""
		A FullIntegration is just two ChimericIntegrations in the same read.
		Assign a ChimericIntegration representing the left junction to self.left
		and a ChimericIntegration representing the right junction to self.right
		"""

		# split viral read into MS ans SM
		vread1, vread2 = self._split_read(self.vread)
		# split host read into SM and MS
		hread1, hread2 = self._split_read(self.hread)
		
		# if host and virus are aligned in same orientation, first part of host read
		# goes with first part of virus read
		if self.hread.is_reverse == self.vread.is_reverse:
				self.left = ChimericIntegration(hread1, vread1, map_thresh =self.map_thresh,
												tol = self.tol, primary = True)
				self.right = ChimericIntegration(hread2, vread2, map_thresh = self.map_thresh, 
												tol = self.tol, primary = True)
		# otherwise, first part of host read goes with second part of virus read and
		# vice versa
		else:
				self.left = ChimericIntegration(hread1, vread2, map_thresh = self.map_thresh,
												tol = self.tol, primary = True)
				self.right = ChimericIntegration(hread2, vread1, map_thresh = self.map_thresh, 
												tol = self.tol, primary = True)						
		
	def get_properties(self):
		"""
		Get a list of properties for each junction.  Properties for each junction is a list
		output by ChimericIntegration.get_properties()
		"""	
		
		props = [self.left.get_properties(), self.right.get_properties()]
		alt_ints = self._get_alt_ints()
		
		# we want to check for rearrangements for the read as a whole, not for individual 
		# junctions
		host_rearrange = self.is_possible_translocation()
		virus_rerrange = self.is_possible_virus_rearrangement()
		
		# set integration properties that are the same for both junctions
		for prop in props:
			prop['Type'] = 'short'
			prop['PossibleVectorRearrangement'] = virus_rerrange
			prop['PossibleHostTranslocation'] = host_rearrange
			
		# assign alternative integrations for each junction
		props[0]['AlternativeInts'] = ";".join([alt_int.left.short_str() for alt_int in alt_ints])
		props[1]['AlternativeInts'] = ";".join([alt_int.right.short_str() for alt_int in alt_ints])			
		
		return props
		
	def _create_segment(self, read, start, end):
		"""
		Create new pysam.AlignedSegment using existing one, using the cigar operations
		between start and end (as in all python, 0-based half-open)
		If there's an insertion at either end of the new CIGAR, convert it to a soft-clip
		
		Get edit distance for this portion of the read by comparing the query and reference
		sequence - getting the query sequence requires the MD tag to be correct, so this 
		needs to be updated if the CIGAR is modified
		"""
		
		assert start >= 0
		assert end > start
		assert end <= len(read.cigartuples)
		
		# must have at least one mapped region in desired part of cigar
		assert 0 in [i[0] for i in read.cigartuples[start:end]]
		
		# get new cigartuples and MD tag
		try:
			md = read.get_tag('MD')
		except KeyError:
			md = ''
			
		# get relevant part of cigartuples (operations between start and end)
		if md != '':
			tmp_cigartuples = AlignmentPool()._combine_MD_with_CIGAR(read.cigartuples, md)
			
			# get just cigartuples
			new_cigartuples = [(op[0], op[1]) for op in tmp_cigartuples[start:end]]
			
			# get MD and edit distance
			md_lst = []
			for op in tmp_cigartuples[start:end]:
				md_lst += op[2]
			
			md = "".join([str(i) for i in md_lst])
		else:
			new_cigartuples = read.cigartuples[start:end]
			
		
		# check there's at least one mapped region in new cigartuples
		assert 0 in [i[0] for i in new_cigartuples]
		
		# if we have insertion at either end, convert to soft clip
		if new_cigartuples[0][0] == 1:
			new_cigartuples[0] = (4, new_cigartuples[0][1])
		if new_cigartuples[-1][0] == 1:
			new_cigartuples[-1] = (4, new_cigartuples[-1][1])
			
		# get an edit distance from md_lst (md tag for this part of the read)
		# don't add CO elements, but instead update these tags for this part of the 
		# read and add to new read
		
		nm = 0
		if md == '':
			# if we don't have an MD tag, just divide the edit distance by the number
			# of mapped regions
			n_mapped = sum([1 for i in read.cigartuples if i[0] == 0])
			nm += math.ceil(read.get_tag('NM') / n_mapped )
		else:
			# check each element in md_lst
			for elem in md_lst:
				try:
					# element is a number (run of matches)
					int(elem)
				except ValueError:
					# element is an insertion
					if elem[0] == '^':
						nm += len(elem) - 1
					# element is a single letter mismatch
					else:
						assert len(elem) == 1
						nm += 1

		# get CO tag for mapped regions between start and stop
		try:
			co = read.get_tag('CO')
		except KeyError:
			co = None
			tmp_co = {}
		
		if co is not None:
			# split co
			tmp_co = AlignmentPool()._split_co(co)
			
			co = ''
			for key, val in tmp_co.items():
				if key >= start and key < end:
					co = co + f"{key-start}:{val[0]}:{val[1]}"
							
		# get subset of query sequence for assigning to new read
		query_start = sum([i[1] for i in read.cigartuples[:start] if i[0] in (0, 1, 4)])
		query_end = sum([i[1] for i in read.cigartuples[:end] if i[0] in (0, 1, 4)])
		new_query = read.query_sequence[query_start:query_end]
#		new_quals = copy.copy(read.query_qualities[query_start:query_end])
		
		# if we're getting a mapped region that isn't the first one, will also need
		# to adjust the mapping position - new_pos

		# get indices of mapped cigar elements
		mapped_idx = [i for i in range(len(read.cigartuples)) if read.cigartuples[i][0] == 0]
		# get start of first mapped block which falls between start and stop
		pos_offset = 0
		for i, map_i in enumerate(mapped_idx):
			# I and D elements combined into mapped regions before the one of interest
			# result in an offset for mapping position of the mapped regions after them
			try:
				pos_offset -= tmp_co[map_i][0] # insertions
				pos_offset += tmp_co[map_i][1] # deletions
			except KeyError:
				pass
			if map_i >= start:
				new_pos = read.get_blocks()[i][0] + pos_offset # don't apply ZP offset now
				assert new_pos >= 0
				assert new_pos <= read.header.get_reference_length(read.reference_name)
				break

		# update tags - keep track of original alignment in OA
		orig_strand = '-' if read.is_reverse else '+'
		OA_update = (read.query_name, str(read.reference_start), orig_strand, 
						read.cigarstring, str(read.mapping_quality), 
						str(read.get_tag('NM')))
		OA_update = ",".join(OA_update) + ";"
		try:
			OA_update = read.get_tag('OA') + OA_update
		except KeyError:
			pass

			
		#create new pysam.AlignedSegment
		a = copy.deepcopy(read)
		
		a.query_sequence = new_query
		a.reference_start = new_pos
		a.cigartuples = new_cigartuples
#		a.query_qualities = new_quals
		
		# update tags
		a.set_tag('NM', nm)
		a.set_tag('OA', OA_update)
		
		# set MD tag if appropriate
		if md != '':
			a.set_tag('MD', md)
			
		if co is not None:
			a.set_tag('CO', co)
		
		return a
	
	def _get_alt_ints(self):
		"""
		Get alternative integrations for the host and viral read - this is different
		to the alternative integrations for each junction
		"""
		
		alt_ints = []
		
		if len(self.hsec) == 0 and len(self.vsec) == 0:
			return alt_ints
		
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
					alt_int = FullIntegration(hst, vrs, map_thresh= self.map_thresh, 
												primary = False)
					alt_ints.append(alt_int)
					
				except AssertionError:
					pass	

		alt_ints = self._swap_self_with_primary(alt_ints, host_alns, virus_alns)
		
		alt_ints = self._filter_alt_ints(alt_ints)

		return alt_ints
	
	def _get_middle_region_read_coords(self, read):
		"""
		Get the part of the read occupied by the middle cigar operation, relative to 
		the original read sequence (i.e. reverse cigar if read is mapped in reverse orientation).
		E.g. for CIGAR 79S76M95S return (79, 155) for forward read and (95, 171) for
		reverse read
		"""
		assert len(read.cigartuples) == 3
		if read.is_reverse:
			start = read.cigartuples[-1][1]
		else:
			start = read.cigartuples[0][1]
			
		stop = start + read.cigartuples[1][1]
		
		assert start < stop
		
		return start, stop
		
	def _get_total_edit_dist(self):
		""" 
		Get total edit distance, which is equal sum of host edit distance, viral edit
		distance, and number of bases in gaps at junctions (if present)
		"""
		
		total = self.left.get_total_edit_dist() + self.right.get_total_edit_dist()
		
		# this double-counts the virus edit distance, so subtract this 
		total -= self.left.get_viral_edit_dist()
		
		return total
		
	def _is_chimeric(self):
		"""
		Check if the host and virus reads look like they contain an integration
		where we can see both junctions.
		
		That is, the host alignment should be mapped/insertion/mapped, and the viral 
		alignment should be clipped/mapped/clipped
		"""	
		
		if None in (self.hread.cigartuples, self.vread.cigartuples):
			return False
		
		# should be mapped/insertion/mapped
		if len(self.hread.cigartuples) < 3:
			return False
			
		# should be clipped/mapped/clipped
		if len(self.vread.cigartuples) < 3:
			return False
		
		# host should be mapped, insertion, mapped			
		host_ops = [i[0] for i in self.hread.cigartuples]
		if host_ops != [0, 1, 0]:
			return False
			
		# virus should be clipped, mapped, clipped
		virus_ops = [i[0] for i in self.vread.cigartuples]
		if virus_ops != [4, 0, 4]:
			return False
		
		# mapped regions must be at least self.map_thresh bp long
		if self.hread.cigartuples[0][1] < self.map_thresh:
			return False
		if self.hread.cigartuples[-1][1] < self.map_thresh:
			return False
		if self.vread.cigartuples[1][1] < self.map_thresh:
			return False
			
		# virus mapped region and host inserted region must overlap
		vstart, vstop = self._get_middle_region_read_coords(self.vread)
		hstart, hstop = self._get_middle_region_read_coords(self.hread)
			
		if vstart > hstop:
			return False
		if hstart > vstop:
			return False
		
		return True
		
	def _split_read(self, read):
		"""
		For a read composed of three cigar operations (SMS or MIM), split into
		two reads with two of the three operations each (SM and MS or MI and IM).
		For reads with insertions, convert these cigar operations to soft-clips
		"""
		assert len(read.cigartuples) == 3
		
		cigarops = [i[0] for i in read.cigartuples]
		assert cigarops == [4,0,4] or cigarops == [0,1,0]
		
		# get first part of read
		first = self._create_segment(read, 0, 2)
		second = self._create_segment(read, 1, 3)

		return first, second
		
	def __str__(self):
		"""
		string representation of an integration
		"""
		
		props = self.get_properties()
		
		props = ["\t".join([str(prop[i]) for i in default_header]) for prop in props]
		
		return "\n".join(props)

def _reverse_complement(seq):
	""" Return reverse complement """
	nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N' , 'X':'X',
			'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n'}
	return ''.join(reversed([nn[i] for i in seq]))	
		
if __name__ == "__main__":
	main()
