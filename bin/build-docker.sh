#!/bin/bash

# Patch dune-geometry
pushd dune/dune-geometry
cp ../../patches/remove_std_call_once.patch .
patch -p1 < remove_std_call_once.patch
popd

# Build using dunecontrol
./dune/dune-common/bin/dunecontrol --opts=./dune/config-docker.opts --module=jupyter-kernel --builddir=$(pwd)/build all
./dune/dune-common/bin/dunecontrol --only=jupyter-kernel --builddir=$(pwd)/build bexec make install_kernelspec
