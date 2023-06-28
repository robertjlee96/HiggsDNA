#!/bin/bash

set -e
pwd

echo "======================================="
echo "Installing HiggsDNA"
echo "======================================="

echo "Fixing dependencies in the image"
conda install -y numba>=0.57.0 llvmlite==0.40.0 numpy>=1.22.0
pip install --upgrade dask-lxplus
cd HiggsDNA
pip3 install .
