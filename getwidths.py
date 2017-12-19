#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# date: 29/apr/2017
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

import os,sys,shutil,subprocess
import ConfigParser
import getopt
import numpy as np
import shapefile as sf
from gdal_utils import *
from osgeo import osr
from scipy.spatial.distance import cdist

import pdb
# pdb.set_trace()

def getwidths(argv):

    opts, args = getopt.getopt(argv,"i:")
    for o, a in opts:
        if o == "-i": inifile  = a

    config = ConfigParser.SafeConfigParser()
    config.read(inifile)

    netf   = str(config.get('getwidths','netf'))
    proj   = str(config.get('getwidths','proj'))
    method = str(config.get('getwidths','method'))
    fwidth = str(config.get('getwidths','fwidth'))
    thresh = np.float64(config.get('getwidths','thresh'))
    output = str(config.get('getwidths','output'))

    fname = output

    w = sf.Writer(sf.POINT)
    w.field('x')
    w.field('y')
    w.field('width')

    # coordinates for bank elevations are based in river network mask
    net   = get_gdal_data(netf)
    geo   = get_gdal_geo(netf)
    iy,ix = np.where(net>0)
    x     = geo[8][ix]
    y     = geo[9][iy]

    for i in range(len(x)):

        print "getwidths.py - " + str(len(x)-i)
        
        xmin  = x[i] - thresh
        ymin  = y[i] - thresh
        xmax  = x[i] + thresh
        ymax  = y[i] + thresh

        width,width_geo = clip_raster(fwidth,xmin,ymin,xmax,ymax)

        if method == 'near':
            mywidth = nearpixel(width,width_geo[8],width_geo[9],np.array([[y[i],x[i]]])) # nearest pixel river

        elif method == 'window_max':
            if width[width>0].size > 0:
                print width.size
                mywidth = width[width>0].max()
            else:
                print "INCREASE THRESHOLD, NO VALUES FOUND IN SEARCHING WINDOW"
                print "PROBLEM AT:" +" "+ "x="+str(x[i]) +" "+ "y="+str(y[i])
                sys.exit()

        # write final value in a shapefile
        w.point(x[i],y[i])
        w.record(x[i],y[i],mywidth)
    
    w.save("%s.shp" % fname)

    # write .prj file
    prj = open("%s.prj" % fname, "w")
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj)
    prj.write(srs.ExportToWkt())
    prj.close()

    fmt    = "GTiff"
    nodata = -9999
    name1  = fname+".shp"
    name2  = fname+".tif"
    subprocess.call(["gdal_rasterize","-a_nodata",str(nodata),"-of",fmt,"-tr",str(geo[6]),str(geo[7]),"-a","width","-a_srs",proj,"-te",str(geo[0]),str(geo[1]),str(geo[2]),str(geo[3]),name1,name2])

def nearpixel(array,ddsx,ddsy,XA):

    """
    Find nearest pixel

    array: array with sourcedata
    ddsx: 1-dim array with longitudes of array
    ddsy: 1-dim array with latitudes of array
    XA: point

    """
    nodata = -9999
    _ds    = np.where(array>0)

    # if there are river pixels in the window
    if _ds[0].size > 0:
        XB  = np.vstack((ddsy[_ds[0]],ddsx[_ds[1]])).T
        ind = np.int(cdist(XA, XB, metric='euclidean').argmin())
        res = array[_ds[0][ind],_ds[1][ind]]
    # otherwise
    else:
        res = nodata

    return res

if __name__ == '__main__':
    getwidths(sys.argv[1:])
