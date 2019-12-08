#!/bin/bash

source ~/.bashrc
conda activate snakemake
export PATH="/scratch1/sco305/conda/bin:$PATH"

snakemake -j 100 -s testAlign.sf --cluster-config cluster.json --rerun-incomplete --use-conda --cluster "sbatch --mem-per-cpu {cluster.mem-per-cpu}  -t {cluster.time} --nodes {cluster.nodes} --cpus-per-task {cluster.cpus-per-task} --output {cluster.output}"