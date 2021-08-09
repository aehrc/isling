# Dockerfile for VSeq-Toolkit

FROM ubuntu:18.04

#  $ docker build . -t szsctt/vseq:latest -t szsctt/vseq:1
#  $ docker run --rm -it --entrypoint /bin/bash szsctt/vseq:latest
#  $ docker push szsctt/vseq:latest
#  $ docker push szsctt/vseq:1

ENV TZ=Australia/Sydney
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
ENV DEBIAN_FRONTEND noninteractive
RUN export DEBIAN_FRONTEND

# Install dependencies and application
RUN apt-get -y update && \
 apt-get -y install apt-utils && \
 apt-get -y install git && \
 apt-get -y install wget gawk && \
 apt-get clean

 # Make work directory
 RUN mkdir -p /var/work &&\
    cd /var/work &&\
    git clone https://github.com/CompMeth/VSeq-Toolkit.git
    
COPY VSeq-Toolkit.txt /var/work/VSeq-Toolkit/
    
ENV VSeqToolkit=/var/work/VSeq-Toolkit

