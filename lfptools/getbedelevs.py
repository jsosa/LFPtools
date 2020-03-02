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
    print('loaded data')

    print('calculating bed from banks and depth')
    bnk['bedelev'] = bnk['elevadj'].astype(np.float32) - dpt['depth'].astype(np.float32)
    print(bnk.keys())
#    print(bnk['bedelev'])

    bed = bnk[['x', 'y', 'geometry','bedelev']]
#    bnk.columns = ['x', 'y', 'bedelev']
    print('Writing out data')
#   Just write bed dataframe to file, rather than creating a new dataframe with the same data
#   mybed = gpd.GeoDataFrame(bnk, crs={'init': 'epsg:4326'}, geometry=[
#                             Point(xy) for xy in zip(bed.x.astype(float), bed.y.astype(float))])
    bed.to_file(output+'.shp')

    nodata = -9999
    fmt = "GTiff"
#    name1 = output
#    name2 = os.path.dirname(output) + '/' + \
#        os.path.basename(output).split('.')[0] + '.tif'
    name1 = output + '.shp'
    name2 = output + '.tif'
    mygeo = gdalutils.get_geo(netf)
    subprocess.call(["gdal_rasterize", "-a_nodata", str(nodata), "-of", fmt,"-ot", "Float32", "-co", "COMPRESS=DEFLATE", "-tr", str(mygeo[6]), str(mygeo[7]), "-a",
                     "bedelev", "-a_srs", proj, "-te", str(mygeo[0]), str(mygeo[1]), str(mygeo[2]), str(mygeo[3]), name1, name2])


if __name__ == '__main__':
    getbedelevs_shell(sys.argv[1:])
