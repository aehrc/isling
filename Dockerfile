
#  $ docker build . -t szsctt/isling:latest -t szsctt/isling:1
#  $ docker run --rm -it szsctt/isling:latest /bin/bash
#  $ docker push szsctt/isling:latest
#  $ docker push szsctt/isling:1

FROM mambaorg/micromamba:0.19.1

USER root

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
ENV TZ=Australia/Sydney
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
ENV DEBIAN_FRONTEND noninteractive
RUN export DEBIAN_FRONTEND

RUN apt-get upate && apt-get install -y libgomp1 && rm -rf /var/lib/apt/lists/*

ADD scripts/consolidate_envs.py /opt/isling/scripts/
ADD envs /opt/isling/envs/
RUN micromamba install -n base -c anaconda pip pyyaml=5.3 -y &&\
	python3 /opt/isling/scripts/consolidate_envs.py /opt/isling/envs/*yml /opt/isling/envs/isling.yml &&\
	micromamba env update -n base -f /opt/isling/envs/isling.yml &&\
	micromamba clean --all -y 	

# include isling scripts, etc
ADD scripts /opt/isling/scripts/
ADD snakemake_rules /opt/isling/snakemake_rules
ADD Snakefile /opt/isling/Snakefile
ADD report_illustrations /opt/isling/report_illustrations

# add test files
ADD test /opt/isling/test
ADD run_tests.sh /opt/isling/

WORKDIR /opt/isling

CMD ./run_tests.sh
