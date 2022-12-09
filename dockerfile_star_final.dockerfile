# Source from : https://github.com/alexdobin/STAR

FROM ubuntu:22.04 as builder

# Upgrade ubuntu, install mandatory dependencies and clean
RUN \
    apt-get update && \
    apt-get upgrade -y -o Dpkg::Options::="--force-confold" && \
    apt-get install -y wget zlib1g-dev libbz2-dev libz-dev g++ gcc make  build-essential autoconf automake gcc perl liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz
RUN tar -xzf 2.7.10b.tar.gz 
RUN cd STAR-2.7.10b/source && make STAR 

RUN mv ./STAR-2.7.10b/source/STAR /usr/bin/.

CMD ["./STAR"]