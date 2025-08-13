# build using:
# docker build . --platform=linux/amd64 -t varunkas/genomics:latest
# docker build . -t varunkas/genomics:latest # build much faster
FROM python:3.12-bookworm

# # required by click
# ENV LC_ALL C.UTF-8
# ENV LANG C.UTF-8

ENV DEBIAN_FRONTEND=noninteractive
ENV PACKAGES="git gcc make g++ libboost-all-dev liblzma-dev libbz2-dev zlib1g-dev libcurl4-openssl-dev curl unzip autoconf apt-transport-https ca-certificates gnupg software-properties-common libncurses5-dev wget"

# update package manager, install essentials, and install python 3
RUN apt-get update && \
  apt-get upgrade -y && \
  apt-get install -y --no-install-recommends ${PACKAGES}

# Install hdf5
# make a new directory
RUN mkdir /tmp/install && \
  cd /tmp/install && \
  wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static && \
  mv bedtools.static /usr/local/bin/bedtools && \
  chmod +x /usr/local/bin/bedtools

RUN mkdir -p /tmp/install && \
  cd /tmp/install && \
  wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph && \
  chmod +x bigWigToBedGraph && \
  mv bigWigToBedGraph /usr/local/bin/

RUN pip install --upgrade pip

WORKDIR /work
# docker push varunkas/genomics:latest
# docker run -it --platform=linux/amd64 varunkas/genomics:latest
