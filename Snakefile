#### pipeline for detecting viral integrations in NGS data ####

# the pipeline detecets viral integrations into a host genome in short reads by 
# identifying chimeric reads, discordant read pairs, and short insertions of virus 
# flanked by host sequence on both sides

# note that only paired-end reads are supported
#
# A yaml file specifies which references are to be used with which datasets
# the pipline lives in the intvi_pipeline folder which should be cloned
# from the git repo

# the steps in the pipeline are as follows:
# 1. preprocessing
#   - convert bam files to fastq, if necessary
#   - either merge R1 and R2 using seqprep, or don't do merging 
#     and just concatenate R1 and R2, as specified in the dsets.yaml config file
# 2. alignments
#   - align all the reads to the virus
#   - remove duplicates
#   - extract only the alinged reads to a fastq
#   - align these reads to the host
# 3. perl scripts
#   - run perl scripts with alignment sam  files as inputs to detect integrations
# 4. postprocessing
#   - apply various types of postprocessing: dedup, filter, mask, annotate
#   - generate ouput xlsx files for each dataset
#   - write output files for visualzation in UCSC browser
#  
# the primary output of the pipeline is the xlsx files for each datasets which
# contain the details of the detected integrations.


#### python modules ####

from glob import glob
from os import path, getcwd
import pandas as pd
import pdb

import sys
sys.path.append(os.path.join(os.getcwd(), "snakemake_rules/"))
from snakemake_rules import make_df
from snakemake_rules import make_reference_dict
from snakemake_rules import make_post_args

# construct dataframe with wildcards and other information about how to run analysis

toDo = make_df(config)

for outdir in set(toDo.loc[:,'outdir']):
	toDo.to_csv(os.path.join(outdir, "analysis_parameters.tsv"), sep = "\t")

# construct dictionary with reference names as keys and reference fastas as values

ref_names = make_reference_dict(toDo)

# construct arguments for postprocess.R script for each dataset

POSTARGS, TOSORT, SORTED = make_post_args(config)


#### global wildcard constraints ####

wildcard_constraints:
	virus = "|".join(set(toDo.loc[:,'virus'])),
	samp = "|".join(set(toDo.loc[:,'sample'])),
	dset = "|".join(set(toDo.loc[:,'dataset'])),
	host = "|".join(set(toDo.loc[:,'host'])),
	align_type = "bwaPaired|bwaSingle",
	outpath = "|".join(set(toDo.loc[:,'outdir']))

#### local rules ####
localrules: all, touch_merged, check_bam_input_is_paired

#### target files ####
summary_files = set()
ucsc_files = set()
merged_bed = set()
for i, row in toDo.iterrows():
	summary_files.add(f"{row['outdir']}/summary/{row['dataset']}.xlsx")
	ucsc_files.add(f"{row['outdir']}/summary/ucsc_bed/{row['dataset']}.post.bed")


rule all:
	input: 
		summary_files,
		ucsc_files,
		expand("{outpath}/{dset}/virus_aligned/{samp}.{virus}.bam",
			zip,
			outpath = toDo.loc[:,'outdir'],
			dset = toDo.loc[:,'dataset'],
			samp = toDo.loc[:,'sample'],
			virus = toDo.loc[:,'virus']
			),
		expand("{outpath}/{dset}/host_aligned/{samp}.{host}.readsFrom{virus}.bam",
			zip,
			outpath = toDo.loc[:,'outdir'],
			dset = toDo.loc[:,'dataset'],
			samp = toDo.loc[:,'sample'],
			virus = toDo.loc[:,'virus'],
			host = toDo.loc[:,'host']
			),
		expand("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.post.merged.bed",
			zip,
			outpath = toDo.loc[:,'outdir'],
			dset = toDo.loc[:,'dataset'],
			samp = toDo.loc[:,'sample'],
			virus = toDo.loc[:,'virus'],
			host = toDo.loc[:,'host']
			)
			

#### read preprocessing ####
include: "snakemake_rules/preprocessing.smk"

#### alignments ####
include: "snakemake_rules/alignment.smk"


#### find integrations ####
include: "snakemake_rules/find_ints.smk"

#### postprocessing ####
include: "snakemake_rules/postprocessing.smk"
