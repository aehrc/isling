
#  $ docker build . -t szsctt/polyidus:latest -t szsctt/polyidus:3
#  $ docker run --rm -it szsctt/polyidus:latest /bin/bash
#  $ docker push szsctt/polyidus:latest
#  $ docker push szsctt/polyidus:3

# cd ~/Documents/Projects/viInt/experiments/expt6_simulations/polyidus
# docker run -v$(pwd):/home --rm szsctt/polyidus:3 python /usr/src/app/src/polyidus.py /home/data/hg38/hg38_bwt2_index /home/data/hpv16/hpv16_bowt_ind --fastq /home/data/fastqfiles/SiHa_R1.fastq.gz /home/data/fastqfiles/SiHa_R2.fastq.gz --outdir /home/polyidusOutput

FROM ubuntu:20.04

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
ENV TZ=Australia/Sydney
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
ENV DEBIAN_FRONTEND noninteractive
RUN export DEBIAN_FRONTEND

RUN apt-get update --fix-missing && apt-get install -y wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git 

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc &&\
    /opt/conda/bin/conda update conda python>3 -y &&\
    /opt/conda/bin/conda clean --all -y 

WORKDIR /usr/src/app

COPY polyidus.yml .

RUN /opt/conda/bin/conda install \
		-c bioconda -c conda-forge \
		bowtie2=2.2 \
		bedtools=2.26 \
		samtools=1.9 \
		pysam=0.15 \
		python=3 \
		psutil=5.4 \
		pandas=0.22 &&\
	/opt/conda/bin/conda clean --all -y 

COPY polyidus .

CMD [ "python3", "src/polyidus.py" ]


