FROM mambaorg/micromamba:0.19.1

USER root

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
ENV TZ=Australia/Sydney
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
ENV DEBIAN_FRONTEND noninteractive
RUN export DEBIAN_FRONTEND

ADD scripts/consolidate_envs.py /opt/isling/scripts/
ADD bwa.yml /tmp/bwa.yml
RUN micromamba install -n base -c anaconda pip pyyaml=5.3 -y &&\
	micromamba env update -n base -f /tmp/bwa.yml &&\
	micromamba clean --all -y 	

