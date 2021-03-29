#
# Copyright (c) 2021, NVIDIA CORPORATION.  All rights reserved.
#
# NVIDIA CORPORATION and its licensors retain all intellectual property
# and proprietary rights in and to this software, related documentation
# and any modifications thereto.  Any use, reproduction, disclosure or
# distribution of this software and related documentation without an express
# license agreement from NVIDIA CORPORATION is strictly prohibited.
#

# Inherit from base NVIDIA CUDA 10.1 image with CUDNN
FROM nvidia/cuda:11.0-cudnn8-runtime-ubuntu18.04

# Set python version to 3.6
RUN apt update && apt install -y \
    python3.6 \
    rsync \
    git \
    python3-pip \
    libz-dev \
    vim \
    nano \
    sudo \
    curl \
    wget \
    tmux \
    htop \
    unzip \
    vcftools

RUN ln -nsf /usr/bin/python3.6 /usr/bin/python

RUN git clone https://github.com/clara-parabricks/gwas-data-simulation-public.git
RUN pip3 install -r gwas-data-simulation-public/requirements.txt

