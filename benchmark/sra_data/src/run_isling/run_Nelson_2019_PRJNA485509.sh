eval "$(conda shell.bash hook)"
conda activate snakemake

module load singularity

cd ../../intvi_pipeline

snakemake --keep-going --rerun-incomplete --use-singularity --profile slurm --jobs 100 --configfile ../config/Nelson_2019_PRJNA485509/isling.yml
