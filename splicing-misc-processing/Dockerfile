# Use ubuntu:latest as parent image
FROM lindsayliang/pipeline_parent

WORKDIR /home

RUN wget http://curl.askapache.com/download/curl-7.48.0.tar.gz \
	&& tar zxvf curl-7.48.0.tar.gz \
	&& cd curl-7.48.0 \
	&& ./configure \
	&& make \
	&& make install \
	&& cd .. \
	&& rm curl-7.48.0.tar.gz


RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
RUN tar xjf htslib-1.9.tar.bz2 \
	&& cd htslib-1.9/ \
	&& ./configure --disable-lzma --enable-plugins --enable-libcurl --prefix=/home/htslib-1.9 \
	&& make  \
	&& make install LDLIBS+=-lcurl LDLIBS+=-lcrypto

RUN wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2
RUN tar xjf samtools-1.3.tar.bz2 \
	&& cd samtools-1.3 \
	&& ./configure --enable-plugins --enable-libcurl --prefix=/home/samtools-1.3 \
	&& make all all-htslib \
	&& make install install-htslib 

RUN pip install HTSeq

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary \
	&& mv bedtools.static.binary bedtools \
	&& chmod a+x bedtools

RUN git clone https://github.com/statgen/libStatGen.git \
	 && cd libStatGen/ \
	 && make all \
	 && cd ../ \
	 && git clone https://github.com/statgen/bamUtil.git \
	 && cd bamUtil/ \
	 && make all \
	 && cd .. 

ENV PATH=$PATH:/home/htslib-1.9/:/home/samtools-1.3:/home/:/home/bamUtil/bin/


COPY . /home
