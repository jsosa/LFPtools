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


def getbedelevs_shell(argv):

    myhelp = '''
LFPtools v0.1

Name
----
getbedelevs

Description
-----------
Get bed elevations by substracting depth from banks

Usage
-----
>> lfp-getbedelevs -i config.txt

Content in config.txt
---------------------
[getbedelevs]
output = Shapefile output file path
netf   = Target mask file path
proj   = Output projection in Proj4 format
bnkf   = Shapefile input bank
dptf   = Shapefile input depth
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

    bnkf = str(config.get('getbedelevs', 'bnkf'))
    dptf = str(config.get('getbedelevs', 'dptf'))
    netf = str(config.get('getbedelevs', 'netf'))
    output = str(config.get('getbedelevs', 'output'))
    proj = str(config.get('getbedelevs', 'proj'))

    getbedelevs(bnkf,dptf,netf,output,proj)

def getbedelevs(bnkf,dptf,netf,output,proj):

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
    getbedelevs_shell(sys.argv[1:])
