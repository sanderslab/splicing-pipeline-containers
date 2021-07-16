FROM lindsayliang/pipeline_parent
WORKDIR /home

ENV DEBIAN_FRONTEND=noninteractive

RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/" >> /etc/apt/sources.list

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    apt-get install apt-transport-https && \
    apt-get update && \
    apt-get install -y r-base r-base-dev

# Install leafcutter prereqs
RUN apt-get update && \
    apt-get install -y \
        gdb \
        libcurl4-openssl-dev \
        libffi-dev \
        liblzma-dev \
        libopenblas-base \
        libopenblas-dev \
        libssh2-1-dev \
        libssl-dev \
        libxml2 \
        libxml2-dev \
        software-properties-common \
        vim


RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
    && tar xjf samtools-1.9.tar.bz2 \
    && rm samtools-1.9.tar.bz2 \
    && cd samtools-1.9/ \
    && ./configure --prefix=/home/samtools_install/ --without-curses --disable-bz2 --disable-lzma \
    && make \
    && make install 

ENV PATH=$PATH:/home/samtools-1.9/

RUN git clone https://github.com/davidaknowles/leafcutter

COPY setup.R /home/
RUN Rscript /home/setup.R
RUN echo "Done setup."

ENV PATH=$PATH:/home/leafcutter/scripts/

RUN pip install sklearn
COPY . /home