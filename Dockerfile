# This Dockerfile tries to be compatible with
# https://mybinder.readthedocs.io/en/latest/tutorials/dockerfile.html#preparing-your-dockerfile
FROM jupyter/base-notebook:584f43f06586

# The gnuplot tool installed with Conda below requires libgl1-mesa-glx.
# Apparently, this cannot be fixed from within Anaconda for horrible reasons:
# https://github.com/conda-forge/pygridgen-feedstock/issues/10
# For all the other dependencies we try to stick with Anaconda for
# consistency and minimal disk space consumption!
USER root
RUN apt update && \
    apt install --no-install-recommends --yes \
      libgl1-mesa-glx && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*
USER ${NB_USER}

# Install conda and its dependencies
RUN conda install -c conda-forge \
        cmake \
        cxxopts \
        fortran-compiler \
        gcc_linux-64 \
        gnuplot \
        make \
        suitesparse \
        superlu \
        xeus-cling && \
    conda clean -a -q -y

# Copying the repository into the Docker container
COPY --chown=${NB_UID} . ${HOME}

# Build the Dune kernel
RUN ./bin/build-docker.sh

# We enable JupyerLab!
ENV JUPYTER_ENABLE_LAB=yes
