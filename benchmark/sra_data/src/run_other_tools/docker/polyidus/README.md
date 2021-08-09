# Polyidus container

A container for [Polyidus](https://github.com/hoffmangroup/polyidus), a tool for identifying viral integration


## Versions

### szsctt/polyidus:1

Contains the requirements specified on the GitHub page bowtie2 (v2.2.6), samtools (v1.9), pysam (v0.8.4), bedtools (v2.26.0), pandas (v0.22.0), numpy (1.18.1), and psutil (v5.4.3), and polyidus code from [my fork of the github repo](https://github.com/suzannesctt/polyidus).

### szsctt/polyidus:2

Contains the updated requirements specified on the GitHub page bowtie2 (v2.2.6), samtools (v1.9), pysam (v0.15.0), bedtools (v2.26.0), pandas (v0.22.0), numpy (1.18.1), and psutil (v5.4.3), and polyidus code from [my fork of the github repo](https://github.com/suzannesctt/polyidus), after updating with the developers fix for outputing read information in the exactHpvIntegrations.tsv file [commit](https://github.com/szsctt/polyidus/commit/bc9ab503ed9f2e8c753b7a2e1c7027977d92ea24)

### szsctt/polyidus:3

I noticed that running version 2 of the container resulted in slightly different results each run.  In contrast, supplying dependencies with either conda or a modulefile that Ondrej Hilinka wrote, results in the same outputs.  Was trying to figure out if it's an issue with the container, or the versions of the dependencies that's the issue.

So made this version with slightly different versions of dependencies - now seems to be reproducible


