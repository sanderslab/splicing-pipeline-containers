FROM lindsayliang/pipeline_parent

RUN git clone https://github.com/chaolinzhanglab/olego.git \
	&& cd olego \
	&& make

ENV PATH=$PATH:/home/olego

COPY . /home