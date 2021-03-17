#!/bin/bash

python -m pip install pytest nbval
cd notebooks
py.test --nbval
