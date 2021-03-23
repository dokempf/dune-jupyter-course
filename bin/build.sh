#!/bin/bash

# Patch dune-geometry
pushd dune/dune-geometry
if [ ! -f "remove_std_call_once.patch" ]; then
    cp ../../patches/remove_std_call_once.patch .
    git apply remove_std_call_once.patch
fi
popd

# Build using dunecontrol
./dune/dune-common/bin/dunecontrol --opts=./dune/config.opts --module=jupyter-kernel --builddir=$(pwd)/build all
./dune/dune-common/bin/dunecontrol --only=jupyter-kernel --builddir=$(pwd)/build bexec make install_kernelspec
