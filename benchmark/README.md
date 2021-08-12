## Benchmarking isling

These folders contain the code to reproduce the results from the isling paper.

The main requirements are `conda` and `singularity`.

You will also need to download the references required to run ViFi from [this google drive link](https://drive.google.com/drive/folders/0ByYcg0axX7udeGFNVWtaUmxrOFk).  These should be placed in this directory (`isling/benchmarking`) - make sure to click 'download all', and then unzip and untar all. This should result in a folder called `Release`, with three folders called `GRCh37`, `GRCh38` and `hg19`.  If your filenames or paths are different, you may need to adjsut the config files accordingly.

### Simulated data

Use the runscript `runme.sh` in the `simulated_data` folder. This script must be run from inside this folder.

#### Downloading references

This step uses Eutils to download some references - the runscript automatically creates a `conda` environment called `sim_isling` that contains these.

Make sure that you have also (manually) downloaded the necessary references for ViFi from the link above.

#### Accuracy - simulation and analysis

The number of cores to use for simulating and analysing data can be set in the `runme.sh` script.  An R script (run in a singularity container).  An Rmarkdown report is run inside a singularity container to generate the figures and tables.

#### Runtime

Again, set the number of cores to use in the `runme.sh` script.  Use as many as possible here (for the data in the manuscript, we used 20). Analysis is again in the form of a Rmarkdown report, run inside a singularity container, which generates the figures and tables.

### SRA data

Use the runscript `runme.sh` in the `sra_data` folder.  This script must be run from inside this folder.

#### Downloading data

Download of SRA data requires the [SRA toolkit](https://github.com/ncbi/sra-tools) and `snakemake`.  The runscript automatically creates a `conda` environment called `snakemake_sra` that contains `sra-tools` and `snakemake`.  Note that the workflow is configured not to use the SRA cache (which you can find the location of using `vdb-config -i`, but rather downloads `.sra` files to a temporary cache in the `data/reads/<dataset>/cache` directory, and then deletes them after creation of `fastq` files.  It is also possible to perform the download step on a different machine with access to a shared filesystem (for example, a cluster where internet access is restricted on some nodes).  To do this, specify the hostname for this machine in the file `config/download/sra.yml`.  Other `snakemake` arguments can be adjusted in the script `src/download_sra/run_get_SRA.sh`.

Make sure that you have also (manually) downloaded the necessary references for ViFi from the link above.

#### Running isling

Snakemake arguments can be adjusted in the script `src/run_isling/run_isling.sh` if necessary.  By default the script runs isling locally, using singularity and 1 core.

#### Running other tools

The other tools for comparison (`ViFi`, `seeksv`, `Polyidus`, `VSeq-Toolkit`) are run using a separate snakemake workflow.

#### Generating figures

An `R` script, run inside a `singularity` container, generates the figures and tables from this part of the manuscript.
