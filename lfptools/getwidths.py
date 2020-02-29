#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: j.sosa@bristol.ac.uk / sosa.jeison@gmail.com

import sys
import subprocess
import configparser
import getopt
import numpy as np
import pandas as pd
import geopandas as gpd
import gdalutils
from lfptools import shapefile
from lfptools import misc_utils
from osgeo import osr


def getwidths_shell(argv):

    myhelp = '''
LFPtools v0.1

Name
----
getwidths

Description
-----------
Retrieve river widths from a data set

Usage
-----
>> lfp-getwidths -i config.txt

Content in config.txt
---------------------
[getwidths]
thresh = Searching window threshold in same units as input data set
output = Shapefile output file path
recf   = `Rec` file path
netf   = Target mask file path
proj   = Output projection in Proj4 format
fwidth = Source width file path GDAL format
method = [const_thresh|var_thresh]
fbankfullq  = Source bankfullq shapefile (Optional, to determine variable threshold)
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

    recf = str(config.get('getwidths', 'recf'))
    netf = str(config.get('getwidths', 'netf'))
    proj = str(config.get('getwidths', 'proj'))
    fwidth = str(config.get('getwidths', 'fwidth'))
    output = str(config.get('getwidths', 'output'))
    thresh = np.float64(config.get('getwidths', 'thresh',-1))
    method = str(config.get('getwidths', 'method','const_thresh'))
    fbankfullq = str(config.get('getwidths', 'fbankfullq',''))

####################################################################    # If there are more NaN than real values, all values in link are equal to 30
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

####################################################################
#
def getwidths(recf,netf, proj, fwidth, output,thresh=-1,method = 'const_thresh',fbankfullq=''):
    if method == 'const_thresh':
        print("    running getwidths.py... constant threshold version")
        getwidths_constthresh(recf, netf, proj, fwidth, output, thresh)
    elif method == 'var_thresh':
        print("    running getwidths.py... variable threshold version")
        # use variable threshold, based on fbankfullq (bankfull q)
        # E.g. use larger search distance for major rivers
        # Could alternatively use accumulation, or strahler order
        getwidths_varthresh(recf,netf, proj, fwidth, output,fbankfullq)

####################################################################
#
def getwidths_varthresh(recf,netf, proj, fwidth, output, fbankfullq):

	# Reading XXX_net.tif file
    geo1 = gdalutils.get_geo(netf)

    bankfullq = gpd.read_file(fbankfullq)
	# bankfullq has name: 'bankfullq'

    # Reading XXX_rec.csv file
    rec = pd.read_csv(recf)
    print('loaded data')

	# x and y resolution (degrees)
    xres = geo1[6]
    yres = geo1[7]
    print('data res',xres,yres)

    width = []
    width = np.ones([len(bankfullq)],dtype=np.float32)*30. # 30 is default value
    for row in bankfullq.itertuples():
        #print(row[0],row[1],row[2],row[3],row[4])
        i = row[0]
        x = float(row[1])
        y = float(row[2])
        bfq = max(float(row[3]),1.)
        # Choose some threshold based on bankfull q (bfq)
        thresh = np.log(bfq)/1000. + bfq/1000000. + 2*abs(xres) + 2*abs(yres)

        # come up with minimum width to search for, based on bankfullq
        # This is designed to prevent assigning
        #width values from the tributaries to the major river channels
        minwidth = bfq/100. + 30


        # Get nearest width from datasource
        # Uses Euclidean distance to find nearest point in source
        # `try` included since it may happen that the width database doesn't
        # contains data in the basin if that is the case all values are assigned
        # a 30 m width


        xmin = x - thresh
        ymin = y - thresh
        xmax = x + thresh
        ymax = y + thresh

        dat, geo = gdalutils.clip_raster(fwidth, xmin, ymin, xmax, ymax)
        try:
            iy, ix = np.where(dat > 30)
        except:
            print('Error: point',i,x,y)
            print('Vals:',bfq,thresh,dat)
            continue
        xdat = geo[8][ix]
        ydat = geo[9][iy]

        try:
            dis, ind = misc_utils.near_euc(xdat, ydat, (x, y))
            val = dat[iy[ind], ix[ind]]
            #width.append(val)
            width[i] = val
        except ValueError:
            #width.append(30.)
            continue

	# Add widths to dataframe, then copy to new dataframe
    #bankfullq['width'] = width
    #widths = bankfullq[['x', 'y', 'geometry','width']]

    rec['width'] = width
    #################################################################
    # Group river network per link
    rec.loc[:, 'width'] = rec.groupby('link').width.apply(check_width)

    # Write out files
    print('Writing out data')
    name1 = output + '.shp'
    #widths.to_file(name1)

    w = shapefile.Writer(shapefile.POINT)
    w.field('x')
    w.field('y')
    w.field('width')
    # Writing .shp resulting file
    for x, y, width in zip(rec['lon'], rec['lat'], rec['width']):
        w.point(x, y)
        w.record(x, y, width)
    w.save("%s.shp" % output)

    # write .prj file
    prj = open("%s.prj" % output, "w")
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj)
    prj.write(srs.ExportToWkt())
    prj.close()

    nodata = -9999
    fmt = "GTiff"
#    name1 = output
#    name2 = os.path.dirname(output) + '/' + \
#        os.path.basename(output).split('.')[0] + '.tif'
    name2 = output + '.tif'
    subprocess.call(["gdal_rasterize", "-a_nodata", str(nodata), "-of", fmt,"-ot", "Float32", "-co", "COMPRESS=DEFLATE", "-tr", str(geo1[6]), str(geo1[7]), "-a",
                     "width", "-a_srs", proj, "-te", str(geo1[0]), str(geo1[1]), str(geo1[2]), str(geo1[3]), name1, name2])





####################################################################
#
def getwidths_constthresh(recf, netf, proj, fwidth, output, thresh):

    print("    running getwidths.py...")
    w = shapefile.Writer(shapefile.POINT)
    w.field('x')
    w.field('y')
    w.field('width')

    # Reading XXX_rec.csv file
    rec = pd.read_csv(recf)

    # Get nearest width from datasource
    # Uses Euclidean distance to find nearest point in source
    # `try` included since it may happen that the width database doesn't
    # contains data in the basin if that is the case all values are assigned
    # a 30 m width
    width = []

    for x, y in zip(rec['lon'], rec['lat']):

        xmin = x - thresh
        ymin = y - thresh
        xmax = x + thresh
        ymax = y + thresh

        dat, geo = gdalutils.clip_raster(fwidth, xmin, ymin, xmax, ymax)
        iy, ix = np.where(dat > 30)
        xdat = geo[8][ix]
        ydat = geo[9][iy]

        try:
            dis, ind = misc_utils.near_euc(xdat, ydat, (x, y))
            val = dat[iy[ind], ix[ind]]
            width.append(val)
        except ValueError:
            width.append(np.nan)

    rec['width'] = width


	#################################################################
    # Group river network per link
    rec.loc[:, 'width'] = rec.groupby('link').width.apply(check_width)

   # Writing .shp resulting file
    for x, y, width in zip(rec['lon'], rec['lat'], rec['width']):
        w.point(x, y)
        w.record(x, y, width)
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
                     "-a", "width", "-a_srs", proj, "-te", str(geo[0]), str(geo[1]), str(geo[2]), str(geo[3]), name1, name2])


if __name__ == '__main__':
    getwidths_shell(sys.argv[1:])
