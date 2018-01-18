#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# date: 29/apr/2017
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

import sys
import subprocess
import ConfigParser
import getopt
import numpy as np
import pandas as pd
import gdal_utils
import misc_utils
import shapefile as sf
from osgeo import osr


def getwidths(argv):

    opts, args = getopt.getopt(argv,"i:")
    for o, a in opts:
        if o == "-i": inifile  = a

    config = ConfigParser.SafeConfigParser()
    config.read(inifile)

    recf   = str(config.get('getslopes','recf'))
    netf   = str(config.get('getwidths','netf'))
    proj   = str(config.get('getwidths','proj'))
    fwidth = str(config.get('getwidths','fwidth'))
    output = str(config.get('getwidths','output'))
    thresh = np.float64(config.get('getwidths','thresh'))

    print "    running getwidths.py..."

    w = sf.Writer(sf.POINT)
    w.field('x')
    w.field('y')
    w.field('width')

    # Reading XXX_rec.csv file
    rec = pd.read_csv(recf)

    # Reading XXX_net.tif file
    dat = gdal_utils.get_gdal_data(fwidth)
    geo = gdal_utils.get_gdal_geo(fwidth)

    iy,ix = np.where(dat>0)
    xdat = geo[8][ix]
    ydat = geo[9][iy]

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

    # Filling NaN, backward and forward
    # Grouping per LINK, then perform operation
    rec.loc[:,'width'] = rec.groupby('link').width.fillna(method='ffill')
    rec.loc[:,'width'] = rec.groupby('link').width.fillna(method='bfill')
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

    geo = gdal_utils.get_gdal_geo(netf)

    fmt    = "GTiff"
    nodata = -9999
    name1  = output+".shp"
    name2  = output+".tif"
    subprocess.call(["gdal_rasterize","-a_nodata",str(nodata),"-of",fmt,"-tr",str(geo[6]),str(geo[7]),"-a","width","-a_srs",proj,"-te",str(geo[0]),str(geo[1]),str(geo[2]),str(geo[3]),name1,name2])

if __name__ == '__main__':
    getwidths(sys.argv[1:])
