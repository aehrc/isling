# Dockerfile for VERSE/VirusFinder2.0

#FROM ubuntu:18.04
FROM ubuntu:14.04


ENV TZ=Australia/Sydney
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
ENV DEBIAN_FRONTEND noninteractive
RUN export DEBIAN_FRONTEND

# Install dependencies and application
RUN apt-get -y update && \
 apt-get -y install apt-utils && \
 apt-get -y install awscli && \
 apt-get -y install environment-modules && \
 apt-get -y install tcl && \
 apt-get -y install git && \
 apt-get -y install perl && \
 apt-get -y install ncbi-blast+ && \
 apt-get -y install bowtie2 && \
 #apt-get -y install samtools && \
 apt-get -y install bwa && \
 #apt-get -y install trinity=2.0.3 && \
 apt-get -y install wget && \
 apt-get -y install make && \
 apt-get -y install gcc mono-mcs && \
 apt-get -y install g++ && \
 apt-get -y install libz-dev && \
 apt-get -y install libncurses5-dev libncursesw5-dev  && \
 apt-get -y install cmake  && \
 apt-get -y install rsync  && \
 apt-get -y install unzip  && \
 #apt-get -y install openjdk-8-jre-headless &&\
 apt-get -y install openjdk-7-jre-headless &&\
 apt-get -y install mailutils &&\
 apt-get clean

 # Make work directory
 RUN mkdir -p /var/work &&\
    cd /var/work

 #Other dependencies
 RUN wget https://vorboss.dl.sourceforge.net/project/samtools/samtools/0.1.18/samtools-0.1.18.tar.bz2 &&\
   tar -xvf samtools-0.1.18.tar.bz2 &&\
   cd samtools-0.1.18 &&\
   make CXXFLAGS=-fPIC CFLAGS=-fPIC CPPFLAGS=-fPIC &&\
   cd .. && rm samtools-0.1.18.tar.bz2

RUN wget https://sourceforge.net/projects/svdetect/files/SVDetect/0.80/SVDetect_r0.8.tar.gz  &&\
     tar -xvf SVDetect_r0.8.tar.gz && \
     rm SVDetect_r0.8.tar.gz

RUN wget https://sourceforge.net/projects/trinityrnaseq/files/PREV_CONTENTS/previous_releases/trinityrnaseq_r2013-02-16.tgz &&\
 	tar -xvf trinityrnaseq_r2013-02-16.tgz &&\
 	cd trinityrnaseq_r2013-02-16 &&\
 	make &&\
 	cd .. && rm trinityrnaseq_r2013-02-16.tgz
 	
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.3.0/bowtie-1.3.0-linux-x86_64.zip &&\
	unzip bowtie-1.3.0-linux-x86_64.zip &&\
	rm bowtie-1.3.0-linux-x86_64.zip

 #RUN wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/v2.11.0/trinityrnaseq-v2.11.0.FULL.tar.gz &&\
 #  tar -xvf trinityrnaseq-v2.11.0.FULL.tar.gz &&\
 #  cd trinityrnaseq-v2.11.0 &&\
 #  make   &&\
 #  make plugins  &&\
 #  make install
 

# Install perl modules
#RUN apt-get install -y cpanminus
RUN apt-get -y install make gcc libexpat1-dev

ENV SAMTOOLS="/samtools-0.1.18/"
RUN cpan local::lib && \
 cpan App::cpanminus && \
 cpanm -n Time::Zone && \
 cpanm -n XML::DOM::XPath && \
 cpanm Bio::SeqIO && \
 cpanm Bio::DB::EUtilities && \
 cpanm Bio::DB::Sam  && \
 cpanm Bio::DB::Sam::Constants && \
 cpanm Config::General &&\
 cpanm Tie::IxHash &&\
 cpanm Parallel::ForkManager


# Export home environments if needed to be done outside config file
ENV JAVA_HOME /usr/lib/jvm/java-7-openjdk-amd64/
RUN export JAVA_HOME
ENV PATH="/samtools-0.1.18:${PATH}"
ENV PATH="/trinityrnaseq_r2013-02-16:${PATH}"
ENV PATH="/bowtie-1.3.0-linux-x86_64:${PATH}"

#Copy tool Code and reads data (CHange path based on local setting)
COPY  VERSE/VirusFinder2.0/ /var/work/VirusFinder2.0/
COPY   data/  /var/work/data/

#Check files in current context
 RUN ls -ltr /trinityrnaseq_r2013-02-16
 RUN which Trinity.pl


 # Create docker run script, this pulls the code from git repo and invokes it
 # Parameters are in Batch environment variables
 #   $FQ1 : Path for reads (e.g. /scratch/kingw/virusFinder/simulation/seq_1.fastq.gz)
 #   $FQ2 : Path for reads(e.g. /scratch/kingw/virusFinder/simulation/seq_2.fastq.gz)



## Create config file
RUN echo 'echo Processing Code..........' > /root/config.txt && \
  #echo 'fastq1 = /var/work/$FQ1' >> /root/config.txt && \
  #echo 'fastq2 = /var/work/$FQ2' >> /root/config.txt && \
  echo 'detect_integration = yes' >> /root/config.txt && \
  echo 'detect_mutation = yes ' >> /root/config.txt && \
  echo 'thread_no = 8 ' >> /root/config.txt && \
  echo 'blastn_bin = /usr/bin/blastn' >> /root/config.txt && \
  echo 'bowtie_bin = /usr/bin/bowtie2' >> /root/config.txt && \
  echo 'bwa_bin = /usr/bin/bwa' >> /root/config.txt && \
  echo 'trinity_script = /trinityrnaseq_r2013-02-16/Trinity.pl' >> /root/config.txt && \
  echo 'SVDetect_dir = /SVDetect_r0.8' >> /root/config.txt && \
  echo 'virus_database = /var/work/data/references/test_AAV.fa' >> /root/config.txt && \
  echo 'bowtie_index_human = /var/work/data/bowtie-index/test_human' >> /root/config.txt && \
  echo 'blastn_index_human = /var/work/data/blast-index/test_human' >> /root/config.txt && \
  echo 'blastn_index_virus = /var/work/data/blast-index/test_AAV' >> /root/config.txt && \
  echo 'min_contig_length = 300' >> /root/config.txt && \
  echo 'blastn_evalue_thrd = 0.05' >> /root/config.txt && \
  echo 'similarity_thrd = 0.8' >> /root/config.txt && \
  echo 'chop_read_length = 25' >> /root/config.txt && \
  echo 'minIdentity = 80' >> /root/config.txt

## Create Entry file
RUN echo 'echo Processing Code..........' > /root/process.sh && \
  echo ' echo fastq1 = /var/work/$FQ1 >> /root/config.txt ' >> /root/process.sh && \
  echo ' echo fastq2 = /var/work/$FQ2 >> /root/config.txt ' >> /root/process.sh && \
  echo ' cat /root/config.txt ' >> /root/process.sh && \
  echo ' chmod -R a+rx /var/work/ ' >> /root/process.sh && \
  echo ' cd /var/work/' >> /root/process.sh && \
  echo ' mkdir results' >> /root/process.sh && \
  echo ' ls -ltr /var/work/VirusFinder2.0/' >> /root/process.sh && \
  echo '/var/work/VirusFinder2.0/VirusFinder.pl -c /root/config.txt ' >> /root/process.sh

#fix permissions
RUN chmod 755 /root/config.txt
RUN chmod 755 /root/process.sh

ENTRYPOINT [ "sh", "-c", "/root/process.sh" ]
