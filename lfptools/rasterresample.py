#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: j.sosa@bristol.ac.uk / sosa.jeison@gmail.com

import os
import sys
import getopt
import configparser
import numpy as np
import multiprocessing as mp
import gdalutils
from osgeo import osr


def rasterresample_shell(argv):

    myhelp = '''
LFPtools v0.1

Name
----
rasterresample

Description
-----------
Resample a DEM by upscaling. It applies a reductions method like
mean, min or meanmin. Outlier detection is also available before running
the reduction method. nproc option defines number of cores to be used when
resampling.

Usage
-----
>> lfp-rasterresample -i config.txt

Content in config.txt
---------------------
[rasterresample]
nproc    = Number of cores to use
outlier  = Outlier detection yes/no
method   = Reduction method mean, min, meanmin
hrnodata = High resolution NODATA value
thresh   = Searching windows threshold
demf     = High resolution DEM
netf     = Target mask file path
output   = Output file, GeoTIFF output
'''

    try:
        opts, args = getopt.getopt(argv, "i:")
        for o, a in opts:
            if o == "-i":
                inifile = a
    except:
        print(myhelp)
        sys.exit(0)

    config = configparser.SafeConfigParser()
    config.read(inifile)

    method = str(config.get('rasterresample', 'method'))
    demf = str(config.get('rasterresample', 'demf'))
    netf = str(config.get('rasterresample', 'netf'))
    output = str(config.get('rasterresample', 'output'))
    outlier = str(config.get('rasterresample', 'outlier'))
    hrnodata = np.float64(config.get('rasterresample', 'hrnodata'))
    thresh = np.float64(config.get('rasterresample', 'thresh'))
    nproc = np.float64(config.get('rasterresample', 'nproc')
                       )  # number of cpus to use

    rasterresample(method,demf,netf,output,outlier,hrnodata,thresh,nproc)

def rasterresample(method,demf,netf,output,outlier,hrnodata,thresh,nproc):

    print("    running rasterresample.py...")

    fname1 = demf
    fname2 = output

    # coordinates for bank elevations are based in river network mask
    net = gdalutils.get_data(netf)
    geo = gdalutils.get_geo(netf)

    # consider all pixels in net30 including river network pixels
    iy, ix = np.where(net > -1)
    x = geo[8][ix]
    y = geo[9][iy]

    # Split x and y in nproc parts
    split_x = np.array_split(x, nproc)
    split_y = np.array_split(y, nproc)

    # Define a queue
    queue = mp.Queue()

    # Setup a list of processes that we want to run
    processes = []
    processes = [mp.Process(target=calc_resampling_mp, args=(
        i, queue, fname1, hrnodata, split_x[i], split_y[i], thresh, outlier, method)) for i in range(len(split_x))]

    # Run processes
    for p in processes:
        p.start()

    # Get process results from the queue
    results = [queue.get() for p in processes]

    # Retrieve results in a particular order
    results.sort()
    results = [r[1] for r in results]

    # Stack results horizontally
    elev = np.hstack(results).reshape(net.shape)

    # Replace NaN by hrnodata
    elev[np.isnan(elev)] = hrnodata

    # elev = calc_resampling(fname1,hrnodata,x,y,ix,iy,thresh,outlier,method)
    gdalutils.write_raster(elev, fname2, geo, "Float32", hrnodata)


def calc_resampling_mp(pos, queue, fname1, hrnodata, x, y, thresh, outlier, method):

    elev = np.ones([len(x)])*hrnodata

    for i in range(len(x)):

        # print("rasterresample.py - " + str(len(x)-i))

        xmin = x[i] - thresh
        ymin = y[i] - thresh
        xmax = x[i] + thresh
        ymax = y[i] + thresh

        dem, dem_geo = gdalutils.clip_raster(fname1, xmin, ymin, xmax, ymax)
        ddem = np.ma.masked_values(dem, hrnodata)
        shape = dem.shape

        # Check for outliers
        if outlier == "yes":
            ddem = check_outlier(dem, ddem, hrnodata, 3.5)

        # Method to applied at every kernel
        if method == "meanmin":
            elev[i] = np.mean([ddem.mean(), ddem.min()])

        elif method == "mean":
            elev[i] = ddem.mean()

        elif method == "min":
            elev[i] = ddem.min()

        else:
            sys.exit('ERROR method not specified')

    queue.put((pos, elev))


def calc_resampling(fname1, hrnodata, x, y, ix, iy, thresh, outlier, method):

    elev = np.ones([len(np.unique(y)), len(np.unique(x))])*hrnodata

    for i in range(len(x)):

        # print("rasterresample.py - " + str(len(x)-i))

        xmin = x[i] - thresh
        ymin = y[i] - thresh
        xmax = x[i] + thresh
        ymax = y[i] + thresh

        dem, dem_geo = gdalutils.clip_raster(fname1, xmin, ymin, xmax, ymax)
        ddem = np.ma.masked_where(dem == hrnodata, dem)
        shape = dem.shape

        if outlier == "yes":
            ddem = check_outlier(dem, ddem, hrnodata, 3.5)
        elev[iy[i], ix[i]] = np.mean([ddem.mean(), ddem.min()])

    return elev


def check_outlier(dem, ddem, hrnodata, thresh):

    shape = dem.shape
    chk = is_outlier(ddem.reshape(-1, 1), thresh)
    arr = np.where(chk == True)
    if arr[0].size > 0:
        dem1 = dem.reshape(-1, 1)
        # dem1_tmp = np.copy(dem1)
        dem1[arr[0]] = hrnodata
        dem2 = dem1.reshape(shape)
        ddem = np.ma.masked_where(dem2 == hrnodata, dem2)

        # # DEBUG DEBUG DEBUG
        # xaxis = np.arange(0,len(dem1)).reshape(-1,1)
        # plt.scatter(xaxis[dem1>hrnodata], dem1[dem1>hrnodata],  color='black')
        # plt.scatter(xaxis[arr[0]], dem1_tmp[arr[0]], color='red')
        # plt.show()

    return ddem


def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """

    if len(points.shape) == 1:
        points = points[:, None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh


if __name__ == '__main__':
    rasterresample_shell(sys.argv[1:])
