## Benchmarking isling

These folders contain the code to reproduce the results from the isling paper.

The main requirements are `conda` and `singularity`.



### Simulated data

### SRA data

Download of SRA data requires the [SRA toolkit](https://github.com/ncbi/sra-tools).  Download a recent release [here](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit), and fill in the path to the `bin` directory in the config files in the file `config/sra.yml`.

If you would like to download the data on a different machine with access to a shared filesystem, specify the hostname for this machine in the file `config/download/sra.yml`




