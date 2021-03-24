# This Dockerfile tries to be compatible with
# https://mybinder.readthedocs.io/en/latest/tutorials/dockerfile.html#preparing-your-dockerfile
FROM jupyter/scipy-notebook:9fe5186aba96

# Install Dune dependencies
USER root
RUN apt update && \
    apt install -yy \
        build-essential \
        cmake \
        git \
        gnuplot \
        libsuitesparse-dev \
        libsuperlu-dev \
        libopenmpi-dev \
        openmpi-bin && \
    rm -rf /var/lib/apt/lists/*

USER ${NB_USER}

# Install conda and its dependencies
RUN conda update --all && \
    conda install xeus-cling cxxopts fortran-compiler -c conda-forge && \
    conda clean -a -q -y

# Copying the repository into the Docker container
COPY --chown=${NB_UID} . ${HOME}

# Build the Dune kernel
RUN ./bin/build.sh

# Make sure that the home directory belongs to the user
USER root
RUN chown -R ${NB_UID} ${HOME}/.local/share/jupyter
USER ${NB_USER}
