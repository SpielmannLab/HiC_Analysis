# build using:
# docker build . --platform=linux/amd64 -t varunkas/fanc:latest
# docker build . -t varunkas/fanc:latest # build much faster
FROM python:3.12-bookworm

# # required by click
# ENV LC_ALL C.UTF-8
# ENV LANG C.UTF-8

ENV DEBIAN_FRONTEND=noninteractive
# ENV PACKAGES="git gcc make g++ libboost-all-dev liblzma-dev libbz2-dev zlib1g-dev libcurl4-openssl-dev curl unzip autoconf apt-transport-https ca-certificates gnupg software-properties-common libncurses5-dev wget"

# # update package manager, install essentials, and install python 3
# RUN apt-get update && \
#   apt-get upgrade -y && \
#   apt-get install -y --no-install-recommends ${PACKAGES}

# RUN mkdir -p /tmp/install/hdf5-build && \
#   cd /tmp/install/hdf5-build
#
# # Install hdf5
# # make a new directory
# RUN wget https://support.hdfgroup.org/releases/hdf5/v1_14/v1_14_6/downloads/hdf5-1.14.6.tar.gz && \
#   tar xzf hdf5-1.*.tar.gz && \
#   rm *.gz
#
# RUN cd hdf5* && \
#   ./configure --prefix="/usr/local" && \
#   make && \
#   make install && \
#   cd / && \
#   rm -rf /tmp/install/hdf5-build

# install fanc latest as per fanc installation instructions - might have to install h5ad before this.
RUN pip install --upgrade pip
RUN pip install "matplotlib==3.7.3" fanc
# RUN pip install "matplotlib==3.7.3" fanc
# RUN pip install numpy==1.26.4
# RUN pip install "fanc==0.9.25"

WORKDIR /work
# docker push varunkas/fanc:latest
