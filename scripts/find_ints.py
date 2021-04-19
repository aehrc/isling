import argparse
import pysam
import pdb

# note python 3.7+ is required

def main():
	
	parser = argparse.ArgumentParser(description='Identify integrations by comparing a host and viral alignment')
	parser.add_argument('--host', help="host alignment bam (query-sorted)", required=True)
	parser.add_argument('--virus', help="virus alignment bam (query-sorted)", required=True)
	parser.add_argument('--verbose', '-v', help="print extra messages?", action='store_true')
	args = parser.parse_args()

	pair = AlignmentFilePair(args.host, args.virus, args.verbose)
	
	pair.find_integrations()
		
	pair.close()
	
		
class AlignmentFilePair:
	""" A class to hold host and viral alignments, and look for integrations/ """
	
	def __init__(self, host, virus, verbose):
		""" Open host and viral alignments, and check for query sorting """
		
		self.host = AlignmentFile(host, verbose)
		self.virus = AlignmentFile(virus, verbose)
		self.verbose = verbose
		
		if self.verbose:
			print(f"using host bam {host} and virus bam {virus}")
			
	def find_integrations(self):
		""" Check for integrations """
		
		# iterate over reads that have alignments to both host and virus
		for alns in self._get_aligns():
			if self.curr_host is None or self.curr_virus is None:
				break
				
			if self.verbose:
				print(f"checking {self.curr_host[0].qname}")
				
			# collect read 1 alignments
			host_r1_primary  = self.curr_host.get_primary_r1()
			virus_r1_primary =  self.curr_virus.get_primary_r1()

			# check if primary read 1 alignments appear to be chimeric	
			if host_r1_primary is not None and virus_r1_primary is not None:
				if self._is_chimeric(host_r1_primary, virus_r1_primary):
					print(f"read 1 chimeric for read {self.curr_host[0].qname}")
					pdb.set_trace()	
					chimera = ChimericRead(host_r1_primary, virus_r1_primary)	
					chimera.get_host_coords()

			# collect read2 alignments
			host_r2_primary = self.curr_host.get_primary_r2()
			virus_r2_primary = self.curr_virus.get_primary_r2()
			
			# check if primary read 1 alignments appear to be chimeric	
			if host_r2_primary is not None and virus_r2_primary is not None:		
				if self._is_chimeric(host_r2_primary, virus_r2_primary):
					print(f"read 2 chimeric for read {self.curr_host[0].qname}")	
				
						
			# if there are alignments that are neither read 1 nor read 2, check if chimeric
			
				
	def _is_chimeric(self, hread, vread):	
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
			ChimericIntegration(hread, vread)
			return True
		except AssertionError:
			return False

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
	
	def get_primary(self):
		""" Return an AlignmentPool with only the primary alignment(s) from a lignmentPool """
		return AlignmentPool([read for read in self if not read.is_supplementary and not read.is_secondary])
		
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

class ChimericIntegration:
	""" A class to store the properties of a simple chimeric integration (ie mapped/clipped and clipped/mapped) """
	
	def __init__(self, host, virus):
		""" Given a host and virus read that are chimeric, calculate the properties of the integration """

		self.hread = host
		self.vread = virus
		
		# simple chimeric read has two CIGAR elements
		assert self._is_chimeric(self.hread, self.vread)
		
		self.chr = host.reference_name
		self.virus = virus.reference_name
		
		self.host_coords = self.get_host_coords()
		
	def get_host_coords(self):
		""" Get coordinates of ambiguous bases in host """
		
		ambig = self.get_ambig_query_coords()
		
		self.hread.get_aligned_pairs(matches_only=False)[ambig[0]]
		pdb.set_trace()
		

	def get_host_edit_dist(self):
		""" Get the edit distance for the host alignment """
		
		return self.hread.get_tag("NM")

	def get_host_seq(self):
		""" Get the part of the read aligned to the host (excluding ambiguous bases) """
		
		return self.hread.query_alignment_sequence
		
	def get_viral_coords(self):
		""" Get coordinates of ambiguous bases in virus """
		
		ambig = self.get_ambig_query_coords()
		if self.gap_or_overlap() == 'overlap':
			pass
			# TODO
		
	def get_viral_edit_dist(self):
		""" Get the edit distance for the viral alignment """
		
		return self.vread.get_tag("NM")
		
	def get_viral_seq(self):
		""" Get the part of the read aligned to the virus (excluding ambiguous bases)"""
		
		return self.vread.query_alignment_sequence
		
	def get_integration_orientation(self):
		""" 
		Get if the junction is host/virus (left side of integration), or virus/host (right side)
		"""
		
		# if host alignment is matched, clipped
		if self.hread.cigartuples[0][0] == 0:
			if self.hread.is_reverse:
				return 'vh'
			else:
				return 'hv'
		
		# if host alignment is clipped, matched
		else:
			if self.hread.is_reverse:
				return 'hv'
			else:
				return 'vh'
		
	def num_ambig_bases(self):
		""" 
		Get the number of ambiguous bases.  This is equal to the number of bases accounted 
		for in both alignments minus the read length.  
		Positive => overlap
		Negative => gap
		Zero => neither
		"""
		
		assert self.hread.infer_read_length() == self.vread.infer_read_length()
		
		total_mapped = self.hread.query_alignment_length + self.vread.query_alignment_length
		
		return total_mapped - self.hread.infer_read_length() 
		
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
			
	def get_ambig_query_coords(self):
		""" Get the coordinates of the ambiguous bases in the read """
		
		# if mapped part of host read is first
		if self.hread.cigartuples[0][0] == 0:
			if self.gap_or_overlap() == 'gap':
				return (
					self.hread.query_alignment_end,
					self.vread.query_alignment_start
				)
			elif self.gap_or_overlap() == 'overlap':
				return (
					self.vread.query_alignment_start,
					self.hread.query_alignment_end
				)
			else:
				return(
					self.hread.query_alignment_end,
					self.hread.query_alignment_end
				)
				
		# if mapped part of virus read is first
		else:
			if self.gap_or_overlap() == 'gap':
				return (
					self.vread.query_alignment_end,
					self.hread.query_alignment_start
				)
			elif self.gap_or_overlap() == 'overlap':
				return (
					self.hread.query_alignment_start,
					self.vread.query_alignment_end
				)
			else:
				return (
					self.hread.query_alignment_start,
					self.hread.query_alignment_start
				)
				
	def _is_chimeric(self, read1, read2):
		"""
		Return True if two alignments are chimeric (occupy complementary parts of the read)
		"""
				
		# if not aligned, can't be chimeric
		if read1.is_unmapped or read2.is_unmapped:
			return False
			
		# viral and host alignment must have only mapped and clipped regions
		if not (len(read2.cigartuples) == 2 and len(read1.cigartuples) == 2):
			return False
		
		# https://pysam.readthedocs.io/en/latest/api.html#api
		# 0 = matched
		# 4 = soft-clipped
		cigarops1 = [i[0] for i in read1.cigartuples]
		cigarops2 = [i[0] for i in read2.cigartuples]
		
		# must have one matched and one soft-clipped part
		if set(cigarops1) != {0, 4}:
			return False
		if set(cigarops2) != {0, 4}:
			return False
			
		# if both forward or both reverse
		if read1.is_reverse == read2.is_reverse:
			# but cigar operations are the same, then no dice
			if cigarops1 == cigarops2:
				return False
		# if aligned in opposite orientations
		else:
			if cigarops1 != cigarops2:
				return False	

		return True
		
		
	def is_possible_translocation(self):
		""" 
		A read is a possible translocation if most of the bases in the read can be accounted
		for by alignments to the host
		"""
		pass
	
	def is_possible_virus_rearrangement(self):
		""" 
		A read is a possible viral rearrangement if most of the bases in the read can be 
		accounted for by alignments to the virus
		"""
		pass
		
		
		
if __name__ == "__main__":
	main()