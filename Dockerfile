ARG FROM_IMAGE=gitlab-registry.cern.ch/batch-team/dask-lxplus/lxdask-cc7:latest
FROM ${FROM_IMAGE}

ARG CLUSTER=lxplus

ADD . . 

RUN echo "=======================================" && \
    echo "Installing HiggsDNA" && \
    echo "=======================================" && \
    if [[ ${CLUSTER} == "lxplus" ]]; then \
    yum -y update && \
    yum -y install git-lfs && \
    echo "Fixing dependencies in the image" && \
    conda install -y numba>=0.57.0 llvmlite==0.40.0 numpy>=1.22.0 && \
    pip install --upgrade dask-lxplus; \
    fi && \
    pip3 install .