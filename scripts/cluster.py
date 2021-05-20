import sys
import argparse
import itertools
import pdb

# cluster integration sites identified with find_ints.py

# find integration sites that cluster together in host and vector.  Integration sites may
# have more than one possible location in the host or vector genomes - integrations
# cluster together if they have at least one possible location that overlaps in both
# host and vector genome

# Essentially, we are making a graph of all integrations. We connect two integrations
# if they have one possible location that overlaps in both host and vector.  We are 
# interested in the connected components in this graph - each is a "Cluster"

# we output the clusters with only one possible location in host and virus to one file,
# and the clusters with multiple possible locations in host or virus to another file

# output coordinates is the range encompassing coordinates common to all ranges included 
# in the cluster.  eg if a cluster consists of the host coordinates 0-1, 0-150 and 0-10,
# the output coordinates will be 0-1.

# the input file should be sorted.  Each integration will therefore only need to be 
# considered alongside it's neighbours in the sorted file, unless it has alternative 
# locations (in which case it may need to be considered alongside all other integrations)

# Two methods for clustering are available - in 'common', output clusters have some common bases in both host and virus genome
# and in 'exact', integration sites are only clustered if they have exactly the same coordinates.
# In both cases, clustered integration sites must have the same orientation (host-virus or virus-host), and 
# viral orientation (integrated in '+' or '-' orientation)

def main(args):
	parser = argparse.ArgumentParser(description = "merge integration sites based on host and virus coordinates")
	parser.add_argument('--input', '-i', help = 'input file (from postprocessing)', required=True)
	parser.add_argument('--output', '-o', help = 'output file', default='merged.txt')
	parser.add_argument('--min-n', '-n', help = 'minimum number of reads to retain a cluster', default=1, type=int)
	parser.add_argument('--cluster-method', '-c', help = 'method for clustering', choices = {'overlap', 'match'}, default='overlap')
	args = parser.parse_args()
	
	print(args)
	
	with IntegrationFile(args.input) as int_file:
		
		int_file.cluster_ints()
			
	

class IntegrationFile(list):
	"""
	A class to read a file containing integrations and find clusters
	"""
	
	def __init__(self, file_name):
		# open file with integrations
		self.file = open(file_name, 'r')
	
	def cluster_ints(self):
		"""
		Read integration file and cluster integrations
		
		Clusters are integrations that overlap/match in both host and vector. 
		"""
		
		self._cluster_ints_host()

			
	
	def _cluster_ints_host(self):
		"""
		First cluster integrations within host
		"""
		
		# first, get header
		try:
			header = next(self.file)
		except StopIteration:
			print(f"integration file is empty!")
			return
			
		self.fields = header.strip().split('\t')
		
		alt_locs = []
		
		# get first integration and add to a cluster
		try:
			line = next(self.file)
		except StopIteration:
			return
		self.append(Cluster(self.fields, line))
		
		# check all integrations in file
		for line in self.file:
			
			# make a new integration object
			integration = Integration(self.fields, line)
			
			# if this integration has only one location, we only need to consider
			# the last cluster
			if integration['AltLocs'] == '':
				# try adding to last cluster, but if it doesn't overlap/match in host
				# then make new cluster
				try:
					self[-1].add_integration(integration)
				except AssertionError:
					self.append(Cluster(self.fields, line))	
			
			# keep integrations with alternative locations for later		
			else:

				alt_locs.append(integration)
			
		# deal with integrations with alternative locations
		print(len(alt_locs))
		pdb.set_trace()
		for integration in alt_locs:
			clusters_map = [i for i in range(len(self)) if self[i]._check_host(integration)]

			print(len(clusters_map))
				
			# if no matches, add new cluster
			if len(clusters_map) == 0:
				self.append(Cluster(self.fields, line))	
			# if one match, add to that cluster
			elif len(clusters_map) == 1:
					elf[clusters_map[0]].add_integration(integration)
			# otherwise, we need to merge the clusters that 	
			else:
					
				pdb.set_trace()
	
	def __enter__(self):
		"""
		Allows the class to be used as a context manager
		"""
		return self
		
	def __exit__(self, type, value, traceback):
		"""
		Allows the class to be used as a context manager
		"""
		self.file.close()


class Cluster(list):
	""" 
	A list of individual integrations that cluster together in host and 
	vector/virus genomes.  Also keeps track of the coordinates spanned by its constituent
	integrations
	"""

	def __init__(self, header, line, type='overlap'):
		"""
		A cluster initially has one integration.

		A cluster consists of integrations that either overlap or have exact matches
		for their coordinates. 
		
		A cluster will have one method for comparing coordinates - this is either 'overlap',
		in which a read is considered for membership in the cluster on the basis of
		overlap with current cluster members in host/viral coordinates, or 'match', 
		in which there must be an exact match between coordinates for all integrations
		in cluster
		
		Initially, a cluster may consist of integrations that overlap/match only in the 
		host.  It may later be split into smaller clusters of integrations that 
		overlap/match in both host and vector.
		"""
		
		assert type in ('overlap', 'match')
		self.type = type
		
		self.append(Integration(header, line))
		
		self.coords = self[0].coords
		
	def add_integration(self, other_int):
		"""
		Add a new integration to the cluster
		"""
		
		# minimally new integration must overlap/match in host
		self._update_coords(other_int)
		
		self.append(other_int)
		
	def _check_host(self, other_int):
		"""
		Check if an Integration object overlaps/matches with the integrations currently in
		this cluster. 
		
		For an integration to overlap/match, it must overlap/match with ALL of the 
		integrations currently in the cluster
		"""
		
		if self.type == 'overlap':
			for self_int in self:
				if not self_int.overlaps_host(other_int):
					return False
				
			return True
		
		else:
			for self_int in self:
				if not self_int.matches_host(other_int):
					return False
				
			return True
		
	def _update_coords(self, other_int):
		"""
		Update self.coords to incorporate all the coordinates for other_int
		"""	
		
		assert self._check_host(other_int)
		
		print([str(i) for i in self.coords])
		
		pdb.set_trace()
		
class Integration(dict):
	"""
	One integration, which may have multiple mapping locations in host and virus genomes
	"""
	
	def __init__(self, header, line):
		"""
		An Integration is initialised from a line from a file containing integrations
		"""
	
		props = line.strip().split("\t")
		assert len(header) == len(props)
		
		for i, field in enumerate(header):
			self[field] = props[i]
			
		# get locations in host and vector
		self.coords = [HostVirusCoordinate(str(self))]		
		
		for alt in self['AltLocs'].split(';'):
			if alt == '':
				continue 
			alt_locs = alt.split(',')
			alt_str = f"{alt_locs[0]},{alt_locs[1]};{alt_locs[2]},{alt_locs[3]};"
			self.coords.append(HostVirusCoordinate(alt_str))
			
			
	def overlaps_host(self, integration):
		"""
		Check if at least one of the host locations for this integration overlaps with at
		least one of the host locations for another integration
		"""
		# check all combinations of coordinates
		for self_coord, other_coord in itertools.product(self.coords, integration.coords):
			# check for overlaps
			if self_coord.overlaps_host(other_coord):
				return True
					
		return False
		
	def matches_host(self, integration):
		"""
		Check if at least one of the host locations for this integration matches with at
		least one of the host locations for another integration
		"""
		# check all combinations of coordinates
		for self_coord, other_coord in itertools.product(self.coords, integration.coords):
			# check for overlaps
			if self_coord.matches_host(other_coord):
				return True
					
		return False	
			
	def _add_coordinates(self, coord_string):
		"""
		Take a string containing host and virus coordinates
		"""
			
	def __str__(self):
		""" String representation is just the coordinates of integration """
		
		host = f"{self['Chr']}:{self['IntStart']}-{self['IntStop']},{self['Orientation']}"
		virus = f"{self['VirusRef']}:{self['VirusStart']}-{self['VirusStop']}"
		virus = virus + "," + self['VirusOrientation']
		
		return f"{host};{virus};{self['ReadID']}"

class HostVirusCoordinate:
	"""
	This class stores a single set of host coordinates and a set of virus coordinates that 
	define a single location for an integration
	"""
	
	__slots__ = 'host', 'virus'
	
	def __init__(self, coord_string):
		"""
		A HostVirusCoordinates object is initialised from a string:
		"chr:start-stop,host_ori;virus:start-stop,virus_ori"
		
		the string may optionally contain a readID which will be ignored:
		"chr:start-stop,host_ori;virus:start-stop,virus_ori;readID"
		"""
		
		coords = coord_string.split(";")
		assert len(coords) == 2 or len(coords) == 3
		
		self.host = Coordinates(coords[0])
		self.virus = Coordinates(coords[1])

	def matches(self, other_coordinate):
		"""
		Check if self coordinates are exactly the same as other_coordinate
		"""
		
		if not self.matches_host(other_coordinate):
			return False
			
		if not self.matches_virus(other_coordinate):
			return False
			
		return True
		
	def matches_host(self, other_coordinate):
		"""
		Check if host coordinates for self match other coordinate host coords exactly
		"""
		
		return self.host.matches(other_coordinate.host)
		
	def matches_virus(self, other_coordinate):
		"""
		Check if virus coordinates for self match other coordinate virus coords exactly
		"""	

		return self.virus.matches(other_coordinate.virus)
		
	def overlaps(self, other_coordinate):
		"""
		Check if both host and virus coordinates overlap for this integration and
		other_coordinate (another HostVirusCoordinate)
		"""
		
		if not self.overlaps_host(other_coordinate):
			return False
			
		if not self.overlaps_virus(other_coordinate):
			return False
			
		return True
		
	def overlaps_host(self, other_coordinate):
		"""
		Check if host coordinates of self overlap with host coordinates of other_coordinate
		"""
		
		return self.host.overlaps(other_coordinate.host)
		
	def overlaps_virus(self, other_coordinate):
		"""
		Check if virus coordinates of self overlap with virus coordinates of other_coordinate
		"""
		
		return self.virus.overlaps(other_coordinate.virus)
		
	def __str__(self):
		"""
		The string representation for this class re-creates the original input string 
		(almost)
		"""
		
		return str(self.host) + ";" + str(self.virus)
		
		
class Coordinates:
	""" 
	This class stores a single location for host and virus, each as a tuple:
	self.ref, self.start, self.stop, self.ori
	
	For consistency between host and virus coordinates objects, convert 'hv' to '+' and
	'vh' ot '-'
	
	There are also methods to check if the coordinates overlap or match exactly
	"""
	
	__slots__ = 'ref', 'start', 'stop', 'ori'
	
	def __init__(self, coord_string):
		"""
		Coordinates are initialised using a string containing host and virus coordinates
		in the format chr:start-stop,ori
		"""
		
		coords = coord_string.split(",")
		
		assert coords[1] in ('hv', 'vh', '+', '-')
		if coords[1] in ('+', '-'):
			self.ori = coords[1]
		elif coords[1] == 'hv':
			self.ori = '+'
		else:
			self.ori = '-'
		
		coords = coords[0].split(":")
		self.ref = coords[0]
		
		self.start, self.stop = [int(i) for i in coords[1].split("-")]
		
		assert self.stop >= self.start
	
	def matches(self, other_coord):
		"""
		Checks if the coordinates in self exactly match the coordinates in other
		"""	
		
		return str(self) == str(other_coord)
		
	def overlaps(self, other_coord):
		"""
		Takes two coordinates tuples (created with self._process_coords) and checks
		if they overlap each other
		"""	
	
		# different chromosome
		if self.ref != other_coord.ref:
			return False
	
		# different orientation
		if self.ori != other_coord.ori:
			return False
			
		# self completely to the left of other_coord
		if self.start > other_coord.stop:
			return False
		
		# self completely to the right of other_coord			
		if other_coord.start > self.stop:
			return False
			
		return True
		 
	def __str__(self):
		"""
		Print all the coordinates for this object
		"""
		
		return f"{self.ref}:{self.start}-{self.stop},{self.ori}"
	
if __name__ == "__main__":
	main(sys.argv)