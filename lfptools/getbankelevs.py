#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: j.sosa@bristol.ac.uk / sosa.jeison@gmail.com

import sys
import getopt
import subprocess
import configparser
import numpy as np
from lfptools import shapefile
import gdalutils
import gdalutils.extras.haversine as haversine
from osgeo import osr
from scipy.ndimage import distance_transform_edt
from scipy.spatial.distance import cdist


def getbankelevs_shell(argv):

    myhelp = '''
LFPtools v0.1

Name
----
getbankelevs

Description
-----------
Get river banks elevations from a high resolution DEM by reductions method like
nearest neighbour, mean, min or meanmin. Additionally, an outlier detection
can be applied to before running the reduction method.

Usage
-----
>> lfp-getbankelevs -i config.txt

Content in config.txt
---------------------
[getbankelevs]
output   = Shapefile output file path
netf     = Target mask file path
proj     = Output projection in Proj4 format
outlier  = Outlier detection yes/no
method   = Reductions available near, mean, min, meanmin
hrnodata = NODATA value for high resolution DEM
thresh   = Saerching threshold in degrees
hrdemf   = High resolution DEM
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

    output = str(config.get('getbankelevs', 'output'))
    netf = str(config.get('getbankelevs', 'netf'))
    hrdemf = str(config.get('getbankelevs', 'hrdemf'))

    try:
        outlier = str(config.get('getbankelevs', 'outlier'))
    except:
        pass

    proj = str(config.get('getbankelevs', 'proj'))
    method = str(config.get('getbankelevs', 'method'))
    hrnodata = np.float64(config.get('getbankelevs', 'hrnodata'))
    thresh = np.float64(config.get('getbankelevs', 'thresh'))

    getbankelevs(output,netf,hrdemf,proj,method,hrnodata,thresh,outlier)

def getbankelevs(output,netf,hrdemf,proj,method,hrnodata,thresh,outlier):

    print("    running getbankelevs.py...")

    fname = output

    w = shapefile.Writer(shapefile.POINT)
    w.field('x')
    w.field('y')
    w.field('elev')

    # coordinates for bank elevations are based in river network mask
    net = gdalutils.get_data(netf)
    geo = gdalutils.get_geo(netf)
    iy, ix = np.where(net > 0)
    x = geo[8][ix]
    y = geo[9][iy]

    for i in range(len(x)):

        xmin = x[i] - thresh
        ymin = y[i] - thresh
        xmax = x[i] + thresh
        ymax = y[i] + thresh

        dem, dem_geo = gdalutils.clip_raster(hrdemf, xmin, ymin, xmax, ymax)
        ddem = np.ma.masked_where(dem == hrnodata, dem)

        if method == 'near':
            nodata = dem_geo[11]
            dfdem = gdalutils.array_to_pandas(dem, dem_geo, nodata, 'gt')
            arr = haversine.haversine_array(np.array(dfdem['y'].values, dtype='float32'), np.float32(
                dfdem['x'].values), np.float32(y[i]), np.float32(x[i]))
            dfdem['dis'] = np.array(arr)
            dfdem.sort_values(by='dis', inplace=True)
            elev = dfdem.iloc[0, 2]

        elif method == 'meanmin':
            if outlier == "yes":
                ddem = check_outlier(dem, ddem, hrnodata, 3.5)
            elev = np.mean([ddem.mean(), ddem.min()])

        elif method == 'mean':
            if outlier == "yes":
                ddem = check_outlier(dem, ddem, hrnodata, 3.5)
            elev = ddem.mean()

        elif method == 'min':
            if outlier == "yes":
                ddem = check_outlier(dem, ddem, hrnodata, 3.5)
            elev = ddem.min()

        # Write final file in a shapefile

        if np.isfinite(elev):
            w.point(x[i], y[i])
            w.record(x[i], y[i], elev)

    w.save("%s.shp" % fname)

    # Write .prj file
    prj = open("%s.prj" % fname, "w")
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj)
    prj.write(srs.ExportToWkt())
    prj.close()

    fmt = "GTiff"
    nodata = -9999
    bnkname1 = output+".shp"
    bnkname2 = output+".tif"
    subprocess.call(["gdal_rasterize", "-a_nodata", str(nodata), "-of", fmt, "-tr", str(geo[6]), str(geo[7]), "-a",
                     "elev", "-a_srs", proj, "-te", str(geo[0]), str(geo[1]), str(geo[2]), str(geo[3]), bnkname1, bnkname2])


def nearivpixel(ddem, rriv, ddsx, ddsy, XA):
    """
    Nearest river pixel when is possible if not
    take value from land
    """
    nodata = -9999
    _ds = np.where(rriv > 0)

    # if there are river pixels in the window
    if _ds[0].size > 0:
        XB = np.vstack((ddsy[_ds[0]], ddsx[_ds[1]])).T
        ind = np.int(cdist(XA, XB, metric='euclidean').argmin())
        elev = ddem[_ds[0][ind], _ds[1][ind]]

    # otherwise take nearest value from land
    elif np.where(rriv == 0)[0].size > 0:
        _ds = np.where(rriv == 0)
        XB = np.vstack((ddsy[_ds[0]], ddsx[_ds[1]])).T
        ind = np.int(cdist(XA, XB, metric='euclidean').argmin())
        elev = ddem[_ds[0][ind], _ds[1][ind]]

    # should be checked
    else:
        elev = nodata

    return elev


def avgrivpixel(ddem, rriv):
    """
    Average the mean and min of river pixels
    """
    nodata = -9999
    _ds = np.where(rriv > 0)

    # if there are river pixels in the window
    if _ds[0].size > 0:
        elev = np.mean([np.ma.masked_where(rriv == 0, ddem).mean(),
                        np.ma.masked_where(rriv == 0, ddem).min()])
    # otherwise
    else:
        elev = nodata

    return elev


def avgedgpixel(ddem, rriv):
    """
    Average the mean and min of edge pixels
    """
    nodata = -9999
    _ds = np.where(rriv > 0)

    # if there are river pixels in the window
    if _ds[0].size > 0:
        euclidis = distance_transform_edt(1-rriv)
        elev = np.mean([np.ma.masked_where(euclidis == 1, ddem).mean(
        ), np.ma.masked_where(euclidis == 1, ddem).min()])
    # otherwise
    else:
        elev = nodata

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
    getbankelevs_shell(sys.argv[1:])
