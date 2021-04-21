import argparse
import pysam
import pdb

# note python 3.7+ is required

def main():
	
	parser = argparse.ArgumentParser(description='Identify integrations by comparing a host and viral alignment')
	parser.add_argument('--host', help="host alignment bam (query-sorted)", required=True)
	parser.add_argument('--virus', help="virus alignment bam (query-sorted)", required=True)
	parser.add_argument('--map-thresh', help="threshold for length of mapped region", default=20, type=int)
	parser.add_argument('--verbose', '-v', help="print extra messages?", action='store_true')
	args = parser.parse_args()

	pair = AlignmentFilePair(args.host, args.virus, args.map_thresh, args.verbose)
	
	pair.find_integrations()
		
	pair.close()
	
		
class AlignmentFilePair:
	""" A class to hold host and viral alignments, and look for integrations """
	
	def __init__(self, host, virus, map_thresh, verbose):
		""" Open host and viral alignments, and check for query sorting """
		
		self.host = AlignmentFile(host, verbose)
		self.virus = AlignmentFile(virus, verbose)
		self.map_thresh = map_thresh
		self.verbose = verbose
		self.ints = []
		
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
			virus_r1_primary =  self.curr_virus.get_primary_r1()

			# check if primary read 1 alignments appear to be chimeric	
			if host_r1_primary is not None and virus_r1_primary is not None:
				integration = self._is_chimeric(host_r1_primary, virus_r1_primary)
				if integration is not None:
					self.ints.append(integration)

			# collect read2 alignments
			host_r2_primary = self.curr_host.get_primary_r2()
			virus_r2_primary = self.curr_virus.get_primary_r2()
			# check if primary read 1 alignments appear to be chimeric	
			if host_r2_primary is not None and virus_r2_primary is not None:
				integration =  self._is_chimeric(host_r2_primary, virus_r2_primary)
				if integration is not None:
					self.ints.append(integration)		

						
			# if there are alignments that are neither read 1 nor read 2, check if chimeric
			host_not12_primary = self.curr_host.get_primary_not_r1_or_r2()
			virus_not12_primary = self.curr_virus.get_primary_not_r1_or_r2()
			if host_not12_primary is not None and virus_not12_primary is not None:
				integration = self._is_chimeric(host_not12_primary, virus_not12_primary)
				if integration is not None:
					self.ints.append(integration)
					
			# check for a discordant pair
			if all([i is not None for i in [host_r1_primary, host_r2_primary, 
											virus_r1_primary, virus_r2_primary]]):
				# check if integration
				integration = self._is_discordant(host_r1_primary, host_r2_primary, 
											virus_r1_primary, virus_r2_primary)
				if integration is not None:
					self.ints.append(integration)
					integration.get_properties()
					
		
		print(f"found {len(self.ints)} integrations")
			
				
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
			return ChimericIntegration(hread, vread, self.map_thresh)
		except AssertionError:
			return None
			
	def _is_discordant(self, hread1, hread2, vread1, vread2):
		"""
		Checks if a pair looks like an discordant integration
		"""
		
		try:
			return DiscordantIntegration(hread1, hread2, vread1, vread2, self.map_thresh)
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
	""" 
	A class to store/calculate the properties of a simple chimeric integration 
	(ie mapped/clipped and clipped/mapped) 
	"""
	
	def __init__(self, host, virus, map_thresh=20):
		""" Given a host and virus read that are chimeric, calculate the properties of the integration """

		self.hread = host
		self.vread = virus
		self.map_thresh = map_thresh
		
		# simple chimeric read has two CIGAR elements
		assert self._is_chimeric(self.hread, self.vread)
		
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
			'ViralAmbiguousLocaiton': self.is_viral_ambig_loc(),
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
		# TODO
		return False
	
	def is_possible_translocation(self):
		""" 
		A read is a possible translocation if most of the bases in the read can be accounted
		for by alignments to the host
		"""
		#TODO
		return False
	
	def is_possible_virus_rearrangement(self):
		""" 
		A read is a possible viral rearrangement if most of the bases in the read can be 
		accounted for by alignments to the virus
		"""
		
		#TODO
		return False
		
	def is_viral_ambig_loc(self):
		"""
		Check if an integration has ambiguous location in viral genome (i.e. viral part of
		read has multiple equivalent alignments to different parts of viral genome)
		"""	

		# TODO
		return False
		
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
					
	def _is_chimeric(self, read1, read2):
		"""
		Return True if two alignments are chimeric (occupy complementary parts of the read)
		Mapped parts of the read must be at least self.map_thresh bp long
		
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
			
		# mapped regions must be at least 20 bp
		if read1.query_alignment_length < self.map_thresh:
			return False
		if read2.query_alignment_length < self.map_thresh:
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

class DiscordantIntegration(ChimericIntegration):
	""" 
	A class to store/calculate the properties of a discordant integration
	(i.e. read 1 mapped to one reference, read2 mapped to a different reference)
	
	This class inherits methods from ChimericIntegration, but re-defines methods that
	need to be implemented differently for a discordant pair
	"""
	
	def __init__(self, host_r1, host_r2, virus_r1, virus_r2, map_thresh=20):
		
		self.hread1 = host_r1
		self.hread2 = host_r2
		self.vread1 = virus_r1
		self.vread2 = virus_r2
		self.map_thresh = map_thresh
		
		# note that self._is_discordant() also assigns self.hread and self.vread
		assert self._is_discordant()
		
		self.get_host_chr()
		self.get_viral_ref()
		
		
		#test = self.get_properties()
		#print(test)
		#pdb.set_trace()
		
	def is_host_ambig_loc(self):
		""" 
		Check if an integration has ambiguous location in host genome (i.e. host part of
		read has multiple equivalent alignments to different parts of host genome)
		"""
		# TODO
		return False
	
	def is_possible_translocation(self):
		""" 
		A read is a possible translocation if most of the bases in the read can be accounted
		for by alignments to the host
		"""
		#TODO
		return False
	
	def is_possible_virus_rearrangement(self):
		""" 
		A read is a possible viral rearrangement if most of the bases in the read can be 
		accounted for by alignments to the virus
		"""
		
		#TODO
		return False
		
	def is_viral_ambig_loc(self):
		"""
		Check if an integration has ambiguous location in viral genome (i.e. viral part of
		read has multiple equivalent alignments to different parts of viral genome)
		"""	

		# TODO
		return False
		

	def gap_or_overlap(self):
		""" For a discordant pair, OverlapType is always discordant """
		return 'discordant'
	
	def get_ambig_seq(self):
		""" Get sequence of ambiguous bases (if present) """
		#TODO
		return None

	def get_host_coords(self):
		""" Get likely coordinates of integration in host chromosome """
		# TODO
		return (None, None)

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
		Get total edit distance, which is equal sum of host edit distance, and viral edit
		distance
		"""
		return self.get_viral_edit_dist() + self.get_host_edit_dist()
	
	def get_viral_coords(self):
		""" Get likely coordinates of integration in viral reference """		
		# TODO
		return (None, None)

	
	def get_viral_orientation(self):
		""" 
		Figure out if the virus is integrated in forward (+) or reverse (-) orientation
		"""
		#TODO
		return None
			
	def get_viral_ref(self):
		""" Get reference name for viral alignment """
		
		v1_map, v2_map = self._is_mapped(self.vread1), self._is_mapped(self.vread2)
		
		assert v1_map != v2_map
		
		if v1_map:
			self.virus = self.vread1.reference_name
		else:
			self.virus = self.vread2.reference_name
		
		
	def _is_discordant(self):
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
		else:
			self.hread = self.hread2
			self.vread = self.vread1
			
		
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
		
if __name__ == "__main__":
	main()