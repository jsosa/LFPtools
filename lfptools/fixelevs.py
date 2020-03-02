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
import pandas as pd
import geopandas as gpd
from lfptools import shapefile
from lfptools import misc_utils
import statsmodels.api as sm
import gdalutils
from osgeo import osr


def fixelevs_shell(argv):

    myhelp = '''
LFPtools v0.1

Name
----
fixelevs

Description
-----------
Fix DEM elevations using either of two methods:
- yamazaki (Yamazaki et al., 2012 J. Hydrol)
- lowless (Locally Weighted Scatterplot Smoothing)

Usage
-----
>> lfp-fixelevs -i config.txt

Content in config.txt
---------------------
[fixelevs]
output = Shapefile output file path
recf   = `Rec` file path
netf   = Target mask file path
proj   = Output projection in Proj4 format
method = yamazaki, lowless
source = Shapefile input file to fix (e.g from lfp-getbankelevs)
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


    """
    This function uses the output from streamnet function from
    TauDEM, specifically the "coord" and "tree" files to adjust
    DEM values from rivers and tributaries for flood using the 
    algorithm bank4flood (1d)

    First create a temporary file where some coordinates points
    have more than 1 value. It happens because when the algorithm is
    applied upstream-downstream and at confluences several values
    are taken for the same coordinate, this error is removed by
    selecting the minimum elevation value.

    """

    source = str(config.get('fixelevs', 'source'))
    output = str(config.get('fixelevs', 'output'))
    netf = str(config.get('fixelevs', 'netf'))
    recf = str(config.get('fixelevs', 'recf'))
    proj = str(config.get('fixelevs', 'proj'))
    method = str(config.get('fixelevs', 'method'))

    fixelevs(source,output,netf,recf,proj,method)

def fixelevs(source,output,netf,recf,proj,method):

    print("    running fixelevs.py...")

    # Reading XXX_net.tif file
    geo = gdalutils.get_geo(netf)

    # Reading XXX_rec.csv file
    rec = pd.read_csv(recf)

    # Reading XXX_bnk.shp file
    bnk_gdf = gpd.read_file(source)

    # Initiate output shapefile
    w = shapefile.Writer(shapefile.POINT)
    w.field('x')
    w.field('y')
    w.field('elevadj')

    # Retrieving bank elevations from XXX_bnk.shp file
    # Values are stored in rec['bnk']
    rec['bnk'] = bnk_gdf['elev'].astype(float)

    # Adjusting bank values, resulting values
    # are stored in rec['bnk_adj']
    # coordinates are grouped by REACH number
    rec['bnk_adj'] = 0
    recgrp = rec.groupby('reach')
    for reach, df in recgrp:
        ids = df.index
        dem = df['bnk']

        # calc bank elevation
        if method == 'yamazaki':
            adjusted_dem = bank4flood(dem)
        elif method == 'lowless':
            adjusted_dem = lowless(dem)
        else:
            sys.exit('Method not recognised')
        rec['bnk_adj'][ids] = adjusted_dem

    # Writing .shp resulting file
    for i in rec.index:
        w.point(rec['lon'][i], rec['lat'][i])
        w.record(rec['lon'][i], rec['lat'][i], rec['bnk_adj'][i])
    w.save("%s.shp" % output)

    # write .prj file
    prj = open("%s.prj" % output, "w")
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj)
    prj.write(srs.ExportToWkt())
    prj.close()

    nodata = -9999
    fmt = "GTiff"
    name1 = output+".shp"
    name2 = output+".tif"
    subprocess.call(["gdal_rasterize", "-a_nodata", str(nodata), "-of", fmt, "-co", "COMPRESS=DEFLATE", "-tr",
                     str(geo[6]), str(geo[7]), "-a", "elevadj", "-a_srs", proj, "-te", str(geo[0]), str(geo[1]), str(geo[2]), str(geo[3]), name1, name2])


def bank4flood(dem):
    """
    Script to adjust river topography following method described
    in Yamazaki et al. (2012, J. Hydrol)

    """

    # TODO:
    # find flat areas larger than 15pixel s in original DEM
    # W=np.ones(1,dem.size)
    # grouped=[(k, sum(1 for i in g)) for k,g in groupby(dem)]
    # ...

    adjusted_dem = np.array(dem)

    for I in range(dem.size-1):  # -1 to avoid error at boundary

        # bug on first and second elevation values
        if adjusted_dem[1] > adjusted_dem[0]:
            adjusted_dem[0] = adjusted_dem[1]

        if adjusted_dem[I+1] > adjusted_dem[I]:

            midind = I
            middem = adjusted_dem[midind]
            vecdem = adjusted_dem[midind:]

            # identify up pixel
            ii = 0
            # look downstream from pixel i and stop when pixel i+1 < midindex
            while vecdem[ii+1] > middem:
                ii = ii+1

                # avoid problems at boundary downstream
                if ii == vecdem.size-1:
                    break

            lastind = midind+ii+1

            zforw = adjusted_dem[midind:lastind]
            zsort = np.sort(zforw)

            zind = []  # indexes
            zmod = []  # adjusted elevation
            lmod = []  # cost function

            for J in range(zsort.size):

                # identify backward pixel
                jj = 1

                # look backward from midindex and stop when pixel i-1>mid-index
                while adjusted_dem[midind-jj] <= zsort[J]:
                    jj = jj+1
                    if jj > midind:
                        break  # avoid problems at boundary upstream

                backind = midind-jj+1

                # extract DEM following backwardindex:forwardindex
                z = adjusted_dem[backind:lastind]

                zind.append(range(backind, lastind))
                # calc adjusted dem for every case
                zmod.append(np.tile(zsort[J], (1, z.size)))
                # calc cost function for every case
                lmod.append(np.sum(np.abs(z-zmod[J])))

            lmin = np.min(lmod)
            imin = np.where(lmod == lmin)[0][0]

            # print("")
            # print("    problem)
            # print("    cost ....." + str(lmod[imin]))
            # print("    indexes .." + str(np.float64(zind[imin])))
            # print("    before ..." + str(adjusted_dem[zind[imin]]))
            # print("    after ...." + str(zmod[imin][0]))
            # print("")

            # final adjusted dem with minimum cost
            adjusted_dem[zind[imin]] = zmod[imin]

    # remove flat banks
    # adjusted_dem2 = avoid_flat_banks(adjusted_dem)

    # # DEBUG
    # fig, ax = plt.subplots()
    # ax.plot(range(dem.size),dem)
    # ax.plot(range(dem.size),adjusted_dem,'--')
    # plt.show()
    # # DEBUG

    return adjusted_dem


def lowless(dem):
    """
    LOWESS (Locally Weighted Scatterplot Smoothing)
    """

    lowess = sm.nonparametric.lowess
    y = dem
    x = np.arange(0, y.size)
    w = lowess(y, x, frac=1/3)
    return w[:, 1]


if __name__ == '__main__':
    fixelevs_shell(sys.argv[1:])
