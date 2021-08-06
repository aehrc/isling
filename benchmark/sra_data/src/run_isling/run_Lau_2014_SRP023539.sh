eval "$(conda shell.bash hook)"
conda activate snakemake

cd ../../

snakemake --keep-going --rerun-incomplete --use-singularity --profile slurm --jobs 100 --configfile benchmarking/sra_data/config/isling/Lau_2014_SRP023539_hg19.yml


