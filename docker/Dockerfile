FROM jsosa/dkr-taudem:latest

MAINTAINER Jeison Sosa <j.sosa@bristol.ac.uk>

RUN conda install -y \
    gdal \
    rtree && \
    pip --no-cache-dir install -U \
    pip \
    cython \
    scipy \
    xarray \
    scikit-learn \
    statsmodels \
    geopandas && \
    conda clean -y --all

RUN pip install \
    git+https://github.com/jsosa/hydroutils.git \
    git+https://github.com/jsosa/gdalutils.git \
    git+https://github.com/jsosa/LFPtools.git \
