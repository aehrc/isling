# Viral integration tools in containers

Put viral integration tools in containers in order to make them easier to use/more portable.  Use with snakemake workflow for tool comparison.

## Test data

The test data is in `data`.  


### Reads
The folder `data/reads` contains reads in fastq format.  These are from simulated integrations of the host into the virus.

### References

The host and viral references (fasta) are in `data/references`.  We're only interested in tools that can use custom references - if they require a particular host or viral reference, don't bother.
The host virus is `data/references/test_human.fa`, and the virus is `data/references/test_AAV.fa`

## Tools

There are many... Listed in order of priority.

[This comparison paper](https://pubmed.ncbi.nlm.nih.gov/30102374/) also contains some information about how they installed VERSE/VirusFinder2 and Seeksv - might be helpful, but not that detailed.

### VERSE/VirusFinder2

[VirusFinder2/VERSE](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-015-0126-6) is a viral integration pipeline that has a ton of dependencies.  The source code can be downloaded [here](https://bioinfo.uth.edu/VirusFinder/).

Was kind of working with conda environment in `VERSE/env`, (it was a while ago that I ran it) but I think from memory there were issues with really big fastq files.   

### Seeksv

[Seeksv](https://academic.oup.com/bioinformatics/article/33/2/184/2525700) looks to more of just a [script](https://github.com/qiukunlong/seeksv).  Probaly just requires samtools picard, bwa, plus the seeksv scripts.

### INSPIIRED

[INSPIIRED](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5363318/) looks like a bunch of R scripts, amongst other thigns.  In theory can be set up using conda, I haven't tried.

[Github](https://github.com/BushmanLab/INSPIIRED)

### ViFi

[ViFi paper](https://pubmed.ncbi.nlm.nih.gov/29579309/)

Mostly done, just need to add my modified scripts in github fork into a new Dockerfile.

Already comes in docker container (no dockerfile provided).  I have a [fork of the github repo](https://github.com/szsctt/ViFi) with some bug fixes.

Uses several environment variables which I need to sort out in the Dockerfile - not done.

### Polyidus

[Polyidus paper](https://www.biorxiv.org/content/10.1101/2020.02.12.942755v2)

Pretty much done.

