#!/bin/bash

source ~/.bashrc
conda activate integration

snakemake -j 100 -s viralIntegrations.sf --cluster-config cluster.json  --rerun-incomplete --cluster "sbatch --mem-per-cpu {cluster.mem-per-cpu}  -t {cluster.time} --mail-user {cluster.mail-user} --nodes {cluster.nodes} --cpus-per-task {cluster.cpus-per-task} --output {cluster.output}"
