# Dockerfile for isling

FROM ubuntu:20.04


#  $ docker build . -t szsctt/isling:latest -t szsctt/isling:1
#  $ docker run --rm -it szsctt/isling:latest /bin/bash
#  $ docker push szsctt/isling:latest
#  $ docker push szsctt/isling:1


ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
ENV TZ=Australia/Sydney
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
ENV DEBIAN_FRONTEND noninteractive
RUN export DEBIAN_FRONTEND

RUN apt-get update --fix-missing && apt-get install -y wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git 

# https://hub.docker.com/r/continuumio/miniconda/dockerfile
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc &&\
    /opt/conda/bin/conda update conda python>3 -y &&\
    /opt/conda/bin/conda clean --all -y 

# install conda stuff
ADD scripts/consolidate_envs.py /opt/isling/scripts/
ADD envs /opt/isling/envs/
RUN /opt/conda/bin/conda install -n base -c anaconda pip pyyaml=5.3 -y &&\
	python3 /opt/isling/scripts/consolidate_envs.py /opt/isling/envs/*yml /opt/isling/envs/isling.yml &&\
	/opt/conda/bin/conda env update -n base -f /opt/isling/envs/isling.yml &&\
	/opt/conda/bin/conda clean --all -y 	

# include isling scripts, etc
ADD scripts /opt/isling/scripts/
ADD snakemake_rules /opt/isling/snakemake_rules
ADD Snakefile /opt/isling/Snakefile

# add test files
ADD test /opt/isling/test
ADD run_tests.sh /opt/isling/

WORKDIR /opt/isling