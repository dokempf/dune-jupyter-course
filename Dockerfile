# This Dockerfile tries to be compatible with
# https://mybinder.readthedocs.io/en/latest/tutorials/dockerfile.html#preparing-your-dockerfile
FROM jupyter/scipy-notebook:9fe5186aba96

# Install Dune dependencies
USER root
RUN apt update && \
    apt upgrade -yy && \
    apt install -yy \
        build-essential \
        cmake \
        git \
        libsuitesparse-dev \
        libsuperlu-dev \
        libopenmpi-dev \
        openmpi-bin
USER ${NB_USER}

# Install conda and its dependencies
RUN wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh && \
    chmod +x Anaconda3-2020.11-Linux-x86_64.sh && \
    ./Anaconda3-2020.11-Linux-x86_64.sh -b && \
    rm Anaconda3-2020.11-Linux-x86_64.sh

RUN conda init
RUN conda update --all
RUN conda install xeus-cling cxxopts -c conda-forge

# Copying the repository into the Docker container
COPY --chown=${NB_UID} . ${HOME}

# Build the Dune kernel
RUN ./bin/build.sh

# Make sure that the home directory belongs to the user
USER root
RUN chown -R ${NB_UID} ${HOME}/.local/share/jupyter
USER ${NB_USER}
