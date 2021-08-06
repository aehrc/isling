eval "$(conda shell.bash hook)"
conda activate snakemake

module load singularity

cd ../../intvi_pipeline

snakemake --keep-going --rerun-incomplete --use-singularity --profile slurm --jobs 100 --configfile ../config/Sung_2012_PRJEB2869/isling_hg19.yml
