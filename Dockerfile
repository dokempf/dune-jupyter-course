# This Dockerfile tries to be compatible with
# https://mybinder.readthedocs.io/en/latest/tutorials/dockerfile.html#preparing-your-dockerfile
FROM jupyter/scipy-notebook:9fe5186aba96

# Install conda and its dependencies
RUN conda update --all && \
    conda install -c conda-forge \
        cmake \
        cxxopts \
        fortran-compiler \
        gcc_linux-64 \
        gnuplot \
        suitesparse \
        superlu \
        xeus-cling && \
    conda clean -a -q -y

# Copying the repository into the Docker container
COPY --chown=${NB_UID} . ${HOME}

# Build the Dune kernel
RUN ./bin/build-docker.sh
