#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# date: 26/jun/2017
# mail: j.sosa@bristol.ac.uk / sosa.jeison@gmail.com

import os
import sys
import getopt
import subprocess
import configparser
import numpy as np
from osgeo import osr
from lfptools import shapefile
import gdalutils
from scipy.spatial.distance import cdist


def getbedelevs(argv):

    opts, args = getopt.getopt(argv,"i:")
    for o, a in opts:
        if o == "-i": inifile  = a

    config = configparser.SafeConfigParser()
    config.read(inifile)
    
    bnkf   = str(config.get('getbedelevs','bnkf'))
    dptf   = str(config.get('getbedelevs','dptf'))
    netf   = str(config.get('getbedelevs','netf'))
    output = str(config.get('getbedelevs','output'))
    proj   = str(config.get('getbedelevs','proj'))

    print("    running getbedelevs.py...")

    nodata = -9999

    fname = output

    w = shapefile.Writer(shapefile.POINT)
    w.field('x')
    w.field('y')
    w.field('bedelev')

    bank  = np.array(shapefile.Reader(bnkf).records(),dtype='float64')
    depth = np.array(shapefile.Reader(dptf).records(),dtype='float64')

    x = bank[:,0]
    y = bank[:,1]

    for i in range(bank.shape[0]):

        dis,ind = near(x,y,np.array([[depth[i,1],depth[i,0]]]))
        if bank[i,2]>nodata and depth[ind,2]>nodata and dis<0.01:
            bed = bank[i,2]-depth[ind,2]
            w.point(x[i],y[i])
            w.record(x[i],y[i],bed)

    # write final value in a shapefile
    w.save("%s.shp" % fname)

    # write .prj file
    prj = open("%s.prj" % fname, "w")
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj)
    prj.write(srs.ExportToWkt())
    prj.close()

    nodata = -9999
    fmt    = "GTiff"
    name1  = fname+".shp"
    name2  = fname+".tif"
    mygeo  = gdalutils.get_geo(netf)
    subprocess.call(["gdal_rasterize","-a_nodata",str(nodata),"-of",fmt,"-tr",str(mygeo[6]),str(mygeo[7]),"-a","bedelev","-a_srs",proj,"-te",str(mygeo[0]),str(mygeo[1]),str(mygeo[2]),str(mygeo[3]),name1,name2])

def near(ddsx,ddsy,XA):

    XB  = np.vstack((ddsy,ddsx)).T
    dis = cdist(XA, XB, metric='euclidean').min()
    ind = cdist(XA, XB, metric='euclidean').argmin()

    return dis,ind

if __name__ == '__main__':
    getbedelevs(sys.argv[1:])
