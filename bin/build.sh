#!/bin/bash

./dune/dune-common/bin/dunecontrol --opts=./dune/config.opts --module=jupyter-kernel --builddir=$(pwd)/build all
