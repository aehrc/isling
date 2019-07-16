#!/bin/bash

source ~/.bashrc
conda activate snakemake

snakemake -j 100 -s rearrange.sf --cluster-config cluster.json  --rerun-incomplete --use-conda --cluster "sbatch --mem-per-cpu {cluster.mem-per-cpu}  -t {cluster.time} --mail-user {cluster.mail-user} --nodes {cluster.nodes} --cpus-per-task {cluster.cpus-per-task} --output {cluster.output}"
