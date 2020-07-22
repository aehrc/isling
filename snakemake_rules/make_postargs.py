from os import path

# construct arguments for postprocess.R script for each dataset
def make_post_args(config):

	POSTARGS = {}
	TOSORT = []
	SORTED = []
	for dataset in config:
		POSTARGS[dataset] = []
		if "post" in config[dataset]:
			for element in config[dataset]["post"]:
				# need to check if this element is a string or a dict
				if type(element) is str:
					# look for keys to be 'filter' or 'dedup'
					if element == "filter":
						POSTARGS[dataset].append("filter")
					elif element == "dedup":
						POSTARGS[dataset].append("dedup")
				# postprocessing types with files specified will be in ordered dictionaries
				elif type(element) is OrderedDict:
					if "mask-exclude" in element.keys():
						for bed in element["mask-exclude"]:
							POSTARGS[dataset].append("mask-exclude")
							POSTARGS[dataset].append(bed)	
					elif "mask-include" in element.keys():
						for bed in element["mask-include"]:
							POSTARGS[dataset].append("mask-include")
							POSTARGS[dataset].append(bed)
					elif "nearest-gtf" in element.keys():
						for gtf in element["nearest-gtf"]:
							sortedgtf = path.splitext(gtf)[0] + ".sorted.gtf"
							POSTARGS[dataset].append("nearest-gtf")
							POSTARGS[dataset].append(sortedgtf)
							if gtf not in TOSORT:
								TOSORT.append(gtf)
								SORTED.append(sortedgtf)
					elif "nearest-bed" in element.keys():
						for bed in element["nearest-bed"]:
							sortedbed = path.splitext(bed)[0] + ".sorted.bed"
							POSTARGS[dataset].append("nearest-bed")
							POSTARGS[dataset].append(sortedbed)
							if bed not in TOSORT:
								TOSORT.append(bed)
								SORTED.append(sortedbed)
					elif "RNA-seq" in element.keys():
						ref = element["genes"]
						sortedref = path.splitext(ref)[0] + ".sorted" + path.splitext(ref)[1]
						if ref not in TOSORT:
							TOSORT.append(ref)
							SORTED.append(sortedref)
						for tsv in element["counts"]:
							POSTARGS[dataset].append("RNA-seq-gtf")
							POSTARGS[dataset].append(sortedref)
							POSTARGS[dataset].append(element["col"])
							POSTARGS[dataset].append(tsv)
		POSTARGS[dataset] = " ".join(POSTARGS[dataset])
		
		return POSTARGS, TOSORT, SORTED

