#!/usr/bin/env python

# inst: university of bristol
# auth: Peter Uhe
# mail: Peter.Uhe@bristol.ac.uk
# Script based on lfptools-getwidths, but run for bankfullq values

import sys
import subprocess
import configparser
import getopt
import numpy as np
import pandas as pd
import gdalutils
from lfptools import shapefile
from lfptools import misc_utils
from osgeo import osr


def getbankfullq_shell(argv):

    myhelp = '''
LFPtools v0.1

Name
----
getbankfullq

Description
-----------
Retrieve bankfull discharge (q) from a data set

Usage
-----
>> lfp-bankfullq -i config.txt

Content in config.txt
---------------------
[getbankfullq]
thresh = Searching window threshold in same units as input data set
output = Shapefile output file path
recf   = `Rec` file path
netf   = Target mask file path
proj   = Output projection in Proj4 format
fbankfullq = Source bankfull Q file path GDAL(TIF) format
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

    recf = str(config.get('getbankfullq', 'recf'))
    netf = str(config.get('getbankfullq', 'netf'))
    proj = str(config.get('getbankfullq', 'proj'))
    fbankfullq = str(config.get('getbankfullq', 'fbankfullq'))
    output = str(config.get('getbankfullq', 'output'))
    thresh = np.float64(config.get('getbankfullq', 'thresh'))

    getbankfullq(recf, netf, proj, fbankfullq, output, thresh)


def getbankfullq(recf, netf, proj, fbankfullq, output, thresh):

    print("    running getbankfullq.py...")

    w = shapefile.Writer(shapefile.POINT)
    w.field('x')
    w.field('y')
    w.field('bankfullq')

    # Reading XXX_rec.csv file
    rec = pd.read_csv(recf)

    # Get nearest bankfullq from datasource
    # Uses Euclidean distance to find nearest point in source
    # `try` included since it may happen that the bankfullq database doesn't
    # contains data in the basin if that is the case all values are assigned
    # 0 Q
    bankfullq = []
    for x, y in zip(rec['lon'], rec['lat']):

        xmin = x - thresh
        ymin = y - thresh
        xmax = x + thresh
        ymax = y + thresh

        dat, geo = gdalutils.clip_raster(fbankfullq, xmin, ymin, xmax, ymax)
        iy, ix = np.where(dat > 0)
        xdat = geo[8][ix]
        ydat = geo[9][iy]

        try:
            dis, ind = misc_utils.near_euc(xdat, ydat, (x, y))
            val = dat[iy[ind], ix[ind]]
            bankfullq.append(val)
        except ValueError:
            bankfullq.append(np.nan)

    rec['bankfullq'] = bankfullq

    # Group river network per link
    # If there are more NaN than real values, all values in link are equal to 30
    # Otherwise, interpolate real values to fill NaNs
    def check_bankfullq(a):
        b = a.copy()
        c = b.isnull()
        falses = c.sum()
        trues = c.count() - falses
        if trues >= falses:
            return a.interpolate(limit_direction='both')
        else:
            b.loc[:] = 0
            return b
    rec.loc[:, 'bankfullq'] = rec.groupby('link').bankfullq.apply(check_bankfullq)

   # Writing .shp resulting file
    for x, y, bankfullq in zip(rec['lon'], rec['lat'], rec['bankfullq']):
        w.point(x, y)
        w.record(x, y, bankfullq)
    w.save("%s.shp" % output)

    # write .prj file
    prj = open("%s.prj" % output, "w")
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj)
    prj.write(srs.ExportToWkt())
    prj.close()

    geo = gdalutils.get_geo(netf)

    fmt = "GTiff"
    nodata = -9999
    name1 = output+".shp"
    name2 = output+".tif"
    subprocess.call(["gdal_rasterize", "-a_nodata", str(nodata), "-of", fmt, "-co", "COMPRESS=DEFLATE","-tr", str(geo[6]), str(geo[7]),
                     "-a", "bankfullq", "-a_srs", proj, "-te", str(geo[0]), str(geo[1]), str(geo[2]), str(geo[3]), name1, name2])


if __name__ == '__main__':
    getbankfullq_shell(sys.argv[1:])
