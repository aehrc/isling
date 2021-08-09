## Benchmarking isling

These folders contain the code to reproduce the results from the isling paper.

The main requirements are `conda` and `singularity`.



### Simulated data

### SRA data

Download of SRA data requires the [SRA toolkit](https://github.com/ncbi/sra-tools) and `snakemake`.  The runscript automatically creates a `conda` environment called `snakemake_sra` that contains `sra-tools` and `snakemake`.  Note that the workflow is configured not to use the SRA cache (which you can find the location of using `vdb-config -i`, but rather downloads `.sra` files to a temporary cache in the `data/reads/<dataset>/cache` directory, and then deletes them after creation of `fastq` files.  It is also possible to perform the download step on a different machine with access to a shared filesystem (for example, a cluster where internet access is restricted on some nodes).  To do this, specify the hostname for this machine in the file `config/download/sra.yml`.  Other `snakemake` arguments can be adjusted in the script `src/download_sra/run_get_SRA.sh`.

### Running isling

Snakemake arguments can be adjusted in the script `src/run_isling/run_isling.sh` if necessary.  By default the script runs isling locally, using singularity and 1 core.

###