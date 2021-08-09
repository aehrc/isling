# Dockerfile for seeksv

FROM ubuntu:20.04

#  $ docker build . -t szsctt/seeksv:latest -t szsctt/seeksv:1
#  $ docker run --rm -it --entrypoint /bin/bash szsctt/seeksv:latest
#  $ docker push szsctt/seeksv:latest
#  $ docker push szsctt/seeksv:1

ENV TZ=Australia/Sydney
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
ENV DEBIAN_FRONTEND noninteractive
RUN export DEBIAN_FRONTEND

# Install dependencies and application
RUN apt-get -y update && \
 apt-get -y install apt-utils && \
 apt-get -y install bwa && \
 apt-get -y install wget && \
 apt-get -y install make && \
 apt-get -y install gcc mono-mcs && \
 apt-get -y install g++ && \
 apt-get -y install libz-dev libbz2-dev liblzma-dev&& \
 apt-get -y install libncurses5-dev libncursesw5-dev  && \
 apt-get -y install cmake  && \
 apt-get -y install rsync  && \
 apt-get -y install unzip  && \
 apt-get -y install default-jre &&\
 apt-get -y install python3 &&\
 apt-get -y install bedtools &&\ 
 apt-get clean

 # Make work directory
 RUN mkdir -p /var/work &&\
    cd /var/work

 #Other dependencies
 RUN wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 &&\
   tar -xvf samtools-1.11.tar.bz2 &&\
   cd samtools-1.11 &&\
   ./configure &&\
   make CXXFLAGS=-fPIC CFLAGS=-fPIC CPPFLAGS=-fPIC &&\
   cd .. && rm samtools-1.11.tar.bz2

RUN mkdir picard &&\
 cd picard &&\
 wget https://github.com/broadinstitute/picard/releases/download/2.23.8/picard.jar

# Export home environments if needed to be done outside config file
ENV JAVA_HOME /usr/bin/java
RUN export JAVA_HOME
ENV PATH="/samtools-1.11:${PATH}"
ENV PICARD /picard/picard.jar
RUN export PICARD

#Check files in current context
 RUN ls -ltr
 RUN which bwa

#Copy tool Code and reads data (CHange path based on local setting)
#COPY  seeksv /var/work/seeksv/
COPY  seeksv/seeksv /var/work/seeksv/
#COPY   data/  /var/work/data/




 # Create docker run script, this pulls the code from git repo and invokes it
 # Parameters are in Batch environment variables
 #   $VAR1 : input bam
 #   $VAR2 : Output prefix


## Create Entry file
RUN echo 'echo Processing Code..........' > /root/process.sh && \
  echo ' chmod -R a+rx /var/work/ ' >> /root/process.sh && \
  echo ' cd /var/work/' >> /root/process.sh && \
  #Example code ran successfully and commented out. Change or update run commands based on variables
  #echo ' /var/work/seeksv/seeksv/seeksv getclip -o normal /var/work/seeksv/example/normal.sort.bam ' >> /root/process.sh && \
  #echo ' /usr/bin/bwa mem /var/work/seeksv/example/reference/example.fa normal.clip.fq.gz | samtools view  -Sb -o normal.clip.bam - ' >> /root/process.sh && \
  #echo ' /var/work/seeksv/seeksv/seeksv getsv normal.clip.bam /var/work/seeksv/example/normal.sort.bam normal.clip.gz normal.sv normal.clipunmap.gz' >> /root/process.sh  && \
  #echo ' ls -ltr /var/work' >> /root/process.sh  && \
  echo ' /var/work/seeksv/seeksv/seeksv getclip -o prefix /var/work/data/$VAR1' >> /root/process.sh && \
  echo ' /usr/bin/bwa mem /var/work/data/references/test_human.fa $VAR2.clip.fq.gz | samtools view  -Sb -o $VAR2.clip.bam - ' >> /root/process.sh && \
  echo ' /var/work/seeksv/seeksv/seeksv getsv $VAR2.clip.bam $VAR1 $VAR2.clip.gz output.sv.txt output.unmapped.clip.fq.gz' >> /root/process.sh


#fix permissions

RUN chmod 755 /root/process.sh

ENTRYPOINT [ "sh", "-c", "/root/process.sh" ]
