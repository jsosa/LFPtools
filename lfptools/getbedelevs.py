#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: j.sosa@bristol.ac.uk / sosa.jeison@gmail.com

import os
import sys
import getopt
import subprocess
import configparser
import numpy as np
from osgeo import osr
import geopandas as gpd
import gdalutils
from shapely.geometry import Point
from scipy.spatial.distance import cdist


def getbedelevs(argv):

    opts, args = getopt.getopt(argv, "i:")
    for o, a in opts:
        if o == "-i":
            inifile = a

    config = configparser.SafeConfigParser()
    config.read(inifile)

    bnkf = str(config.get('getbedelevs', 'bnkf'))
    dptf = str(config.get('getbedelevs', 'dptf'))
    netf = str(config.get('getbedelevs', 'netf'))
    output = str(config.get('getbedelevs', 'output'))
    proj = str(config.get('getbedelevs', 'proj'))

    print("    running getbedelevs.py...")

    bnk = gpd.read_file(bnkf)
    dpt = gpd.read_file(dptf)

    bnk_buf = gpd.GeoDataFrame(
        bnk, crs={'init': 'epsg:4326'}, geometry=bnk.buffer(0.001))

    bed = gpd.sjoin(bnk_buf, dpt, op='intersects')

    bed['bedelev'] = bed['elevadj'].astype(float) - bed['depth'].astype(float)

    bed = bed[['x_left', 'y_left', 'bedelev']]
    bed.columns = ['x', 'y', 'bedelev']

    mybed = gpd.GeoDataFrame(bed, crs={'init': 'epsg:4326'}, geometry=[
                             Point(xy) for xy in zip(bed.x.astype(float), bed.y.astype(float))])

    mybed.to_file(output)

    nodata = -9999
    fmt = "GTiff"
    name1 = output
    name2 = os.path.dirname(output) + '/' + \
        os.path.basename(output).split('.')[0] + '.tif'
    mygeo = gdalutils.get_geo(netf)
    subprocess.call(["gdal_rasterize", "-a_nodata", str(nodata), "-of", fmt, "-tr", str(mygeo[6]), str(mygeo[7]), "-a",
                     "bedelev", "-a_srs", proj, "-te", str(mygeo[0]), str(mygeo[1]), str(mygeo[2]), str(mygeo[3]), name1, name2])


if __name__ == '__main__':
    getbedelevs(sys.argv[1:])
