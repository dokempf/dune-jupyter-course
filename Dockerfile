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
RUN mkdir -p /opt/dune && chown -R ${NB_USER} /opt/dune
USER ${NB_USER}

# Copying the repository into the Docker container
COPY --chown=${NB_UID} . /opt/dune

# Install Conda environment
RUN conda env update -n base --file /opt/dune/environment.yml && \
    conda clean -a -q -y

# Build the Dune kernel
WORKDIR /opt/dune
RUN ./bin/build-docker.sh
WORKDIR ${HOME}

# Populate the home directory with our notebooks.
RUN cp -r /opt/dune/notebooks/* ${HOME}
