import sys
import argparse
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

# the input file need not be sorted, since 

# Two methods for clustering are available - in 'common', output clusters have some common bases in both host and virus genome
# and in 'exact', integration sites are only clustered if they have exactly the same coordinates.
# In both cases, clustered integration sites must have the same orientation (host-virus or virus-host), and 
# viral orientation (integrated in '+' or '-' orientation)

def main(args):
	parser = argparse.ArgumentParser(description = "merge integration sites based on host and virus coordinates")
	parser.add_argument('--input', '-i', help = 'input file (from postprocessing)', required=True)
	parser.add_argument('--output', '-o', help = 'output file', default='merged.txt')
	parser.add_argument('--min-n', '-n', help = 'minimum number of reads to retain a cluster', default=1, type=int)
	parser.add_argument('--cluster-method', '-c', help = 'method for clustering', choices = {'common', 'exact'}, default='common')
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
		
		Clusters are integrations that overlap in both host and vector.  Consider
		integrations that 
		"""
		
		# first, get header
		try:
			header = next(self.file)
		except StopIteration:
			print(f"integration file is empty!")
			return
			
		self.fields = header.strip().split('\t')
		
		for line in self.file:
			
			integration = Integration(self.fields, line)
			self.append(integration)
			
			
	
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
		
	def _check_overlaps(self, line):
		pass


class Cluster(list):
	""" 
	A list of individual integrations that cluster together in host and 
	vector/virus genomes
	"""

	def __init__(self, header, line):
		"""
		A cluster initially has one integration
		"""
		
		self.append(Integration(header, line))
		
	def overlap(self, integration, type='overlap'):
		"""
		Check if an Integration object overlaps with the integrations currently in
		this cluster. There are two types of overlap we consider:
		'overlap': if one of the locations for the integration overlaps in both host
		and vector with one of the locations in the cluster
		"""
		
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
		self.coords = [Coordinates(str(self))]
		
		
		for alt in self['AltLocs'].split(';'):
			if alt == '':
				continue 
			alt_locs = alt.split(',')
			alt_str = f"{alt_locs[0]},{alt_locs[1]};{alt_locs[2]},{alt_locs[3]};"
			self.coords.append(Coordinates(alt_str))
			
	def __str__(self):
		""" String representation is just the coordinates of integration """
		
		host = f"{self['Chr']}:{self['IntStart']}-{self['IntStop']},{self['Orientation']}"
		virus = f"{self['VirusRef']}:{self['VirusStart']}-{self['VirusStop']}"
		virus = virus + "," + self['VirusOrientation']
		
		return f"{host};{virus};{self['ReadID']}"
		
class Coordinates:
	""" A class that stores a set of host and virus coordinates for an integration """
	
	__slots__ = 'host', 'virus'
	
	def __init__(self, coord_string):
		"""
		Coordinates are initialised using a string containing host and virus coordinates
		in the format chr:start-stop,host_ori;virus:start-stop,virus_ori;readID
		"""

		host, virus, _ = coord_string.split(";")
		
		self.host = self._process_coords(host)
		self.virus = self._process_coords(virus)
		
	def matches(self, coordinates):
		"""
		Check if there is an exact matche between self and another Coordiantes object
		"""
		
		if not self._matches_single(self.host, coordinates.host):
			return False
			
		if not self._matches_single(self.virus, coordinates.virus):
			return False
			
		return True

	def overlaps(self, coordinates):
		"""
		Check if there is an overlap between self and another Coordinates object
		"""
		
		if not self._overlaps_single(self.host, coordinates.host):
			return False
			
		if not self._overlaps_single(self.virus, coordinates.virus):
			return False
			
		return True
	
	def _matches_single(self, coord_tuple1, coord_tuple2):
		"""
		Takes two coordinates tuples (created with self._process_coords) and checks
		if they exactly match each other
		"""	
		
		assert len(coord_tuple1) == len(coord_tuple2)
		for i in range(len(coord_tuple1)):
			assert type(coord_tuple1[i]) == type(coord_tuple2[i])
		
		return coord_tuple1 == coord_tuple2
		
	def _overlaps_single(self, coord_tuple1, coord_tuple2):
		"""
		Takes two coordinates tuples (created with self._process_coords) and checks
		if they overlap each other
		"""	
		
		ref1, start1, stop1, ori1 = coord_tuple1
		ref2, start2, stop2, ori2 = coord_tuple2
		
		assert stop1 >= start1
		assert stop2 >= start2

	
		# different chromosome
		if ref1 != ref2:
			return False
	
		# different orientation
		if ori1 != ori2:
			return False
			
		# interval 1 completely to the left of interval 2
		if start1 > stop2:
			return False
			
		if start2 > stop1:
			return False
			
		return True
	
	def _process_coords(self, coord_string):
		"""
		Take a string containing coordaintes for the host or vector in the format 
		chr:start-stop,ori
		Return a tuple containing the chromosome/virus, start, stop and orientation
		"""
		
		coords = coord_string.split(",")
		
		ori = coords[1]
		assert ori in ('hv', 'vh', '+', '-')
		
		coords = coords[0].split(":")
		chr = coords[0]
		
		start, stop = coords[1].split("-")
		
		return (chr, int(start), int(stop), ori)
		 
	
if __name__ == "__main__":
	main(sys.argv)