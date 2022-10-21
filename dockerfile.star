FROM ubuntu:22.04

RUN apt-get update --fix-missing \
&& apt-get install -y wget gcc g++ make libbz2-dev zlib1g zlib1g-dev \
&& cd /usr/local \
&& wget https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz \
&& tar -xzf 2.7.10a.tar.gz \
&& cd STAR-2.7.10a/source && make STAR 

RUN mkdir /data
WORKDIR /data/
RUN mv /usr/local/STAR-2.7.10a/source/STAR /usr/local/bin/.

CMD ["STAR"]