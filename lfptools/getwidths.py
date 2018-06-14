#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# date: 29/apr/2017
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

import sys
import subprocess
import configparser
import getopt
import numpy as np
import pandas as pd
import gdalutils
from lfptools import shapefile
from lfptools import misc_utils
from lfptools.outlier import outlier
from osgeo import osr


def getwidths(argv):

    opts, args = getopt.getopt(argv,"i:")
    for o, a in opts:
        if o == "-i": inifile  = a

    config = configparser.SafeConfigParser()
    config.read(inifile)

    recf   = str(config.get('getwidths','recf'))
    netf   = str(config.get('getwidths','netf'))
    proj   = str(config.get('getwidths','proj'))
    fwidth = str(config.get('getwidths','fwidth'))
    output = str(config.get('getwidths','output'))
    thresh = np.float64(config.get('getwidths','thresh'))

    print("    running getwidths.py...")

    w = shapefile.Writer(shapefile.POINT)
    w.field('x')
    w.field('y')
    w.field('width')

    # Reading XXX_rec.csv file
    rec = pd.read_csv(recf)

    # Reading width source file
    dat   = gdalutils.get_data(fwidth)
    geo   = gdalutils.get_geo(fwidth)
    iy,ix = np.where(dat>30)
    xdat  = geo[8][ix]
    ydat  = geo[9][iy]

    # Get nearest width from datasource
    # Uses Euclidean distance to find nearest point in source
    width = []
    for x,y in zip(rec['lon'],rec['lat']):
        dis,ind = misc_utils.near_euc(xdat,ydat,(x,y))
        if dis <= thresh:
            val = dat[iy[ind],ix[ind]]
            width.append(val)
        else:
            width.append(np.nan)
    rec['width'] = width

    # Group river network per link
    # If there are more NaN than real values, all values in link are equal to 30
    # Otherwise, interpolate real values to fill NaNs
    def check_width(a):
        b = a.copy()
        c = b.isnull()
        falses = c.sum()
        trues = c.count() - falses
        if trues >= falses:
            return a.interpolate(limit_direction='both')
        else:
            b.loc[:] = 30
            return b
    rec.loc[:,'width'] = rec.groupby('link').width.apply(check_width)

    # Saving width data in csv
    rec.to_csv(output+".csv")

   # Writing .shp resulting file
    for x,y,width in zip(rec['lon'],rec['lat'],rec['width']):
        w.point(x,y)
        w.record(x,y,width)
    w.save("%s.shp" % output)

    # write .prj file
    prj = open("%s.prj" % output, "w")
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj)
    prj.write(srs.ExportToWkt())
    prj.close()

    geo = gdalutils.get_geo(netf)

    fmt    = "GTiff"
    nodata = -9999
    name1  = output+".shp"
    name2  = output+".tif"
    subprocess.call(["gdal_rasterize","-a_nodata",str(nodata),"-of",fmt,"-tr",str(geo[6]),str(geo[7]),"-a","width","-a_srs",proj,"-te",str(geo[0]),str(geo[1]),str(geo[2]),str(geo[3]),name1,name2])

    # # Fix dataset
    # shpf = name1
    # labl = 'width'
    # outf = output
    # outlier(recf,proj,netf,shpf,labl,outf,'qua')

if __name__ == '__main__':
    getwidths(sys.argv[1:])
