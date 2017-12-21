#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# date: 22/apr/2017
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

import os
import sys
import pdb # pdb.set_trace()
import getopt
import shutil
import subprocess
import ConfigParser
import numpy as np
import pandas as pd
import shapefile as sf
from osgeo import osr
from gdal_utils import *

def getinflows(argv):

    opts, args = getopt.getopt(argv,"i:")
    for o, a in opts:
        if o == "-i": inifile = a

    config = ConfigParser.SafeConfigParser()
    config.read(inifile)

    method = str(config.get('getinflows','method'))
    treef  = str(config.get('getinflows','treef'))
    netf   = str(config.get('getinflows','netf'))
    coordf = str(config.get('getinflows','coordf'))
    proj   = str(config.get('getinflows','proj'))
    thresh = np.float64(config.get('getinflows','thresh'))
    output = str(config.get('getinflows','output'))

    mygeo = get_gdal_geo(netf)

    fname = output

    w = sf.Writer(sf.POINT)
    w.field('x')
    w.field('y')
    w.field('linkno')
    w.field('inflow')

    # loading coord and tree files
    coord_df = read_coord(coordf)
    tree_df  = read_tree(treef)
    nlinks   = tree_df['linkno'].size

    # creating empty lists to store final lons and lats
    mylon = []
    mylat = []
    mylno = []

    # iterate over links
    for i in range(nlinks):

        print "getinflows.py - " + str(nlinks-i)

        linkno = tree_df.linkno[i]
        start  = tree_df.start[i]
        end    = tree_df.end[i] + 1 # +1 because inlusive slice

        """ Locations in the middle of every reach with no upstream reach"""

        # if tree[i,4] == -1.:
        #     midloc = int(np.median(np.arange(start,end)))
        #     lon = coord[midloc,0]
        #     lat = coord[midloc,1]
        #     mylon.append(lon)
        #     mylat.append(lat)

        """ Locations at the beginning of every reach"""

        # lon = coord[start,0]
        # lat = coord[start,1]
        # mylon.append(lon)
        # mylat.append(lat)

        """ Locations in middle of every reach """

        # midloc = int(np.median(np.arange(start,end)))
        # lon    = coord[midloc,0]
        # lat    = coord[midloc,1]
        # mylon.append(lon)
        # mylat.append(lat)

        """ Locations every 5KM """

        lons = coord_df.lon[start:end].values
        lats = coord_df.lat[start:end].values

        # get coordinates first point
        firstlat = lats[0]
        firstlon = lons[0]

        # iterate along lat, lon vectors
        for j in range(len(lons)-1): # -1, because error outbound

            if method == 'haversine':
                dis = haversine([firstlat,firstlon],[lats[j+1],lons[j+1]])

            # check if point is larger than a threshold
            if dis >= thresh:
                mylno.append(linkno)
                mylon.append(firstlon)
                mylat.append(firstlat)
                firstlat = lats[j+1]
                firstlon = lons[j+1]

    # create a pandas dataframe
    d  = {'x': mylon, 'y':mylat, 'linkno': mylno,}
    df = pd.DataFrame(data=d, index=range(len(mylno)))

    # clean links
    grouped = df.groupby('linkno')

    df_new = pd.DataFrame()
    for name,group in grouped:
        if group['linkno'].size >= 3:
            df_new = pd.concat((df_new,group))

    # final dataframe
    df_new = df_new.reset_index(drop=True)

    # write coordinate points in shapefile
    for i in range(df_new['linkno'].size):
        w.point(df_new['x'][i],df_new['y'][i])
        w.record(df_new['x'][i],df_new['y'][i],df_new['linkno'][i],df_new.index.values[i])
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
    subprocess.call(["gdal_rasterize","-a_nodata",str(nodata),"-of",fmt,"-tr",str(mygeo[6]),str(mygeo[7]),"-a","inflow","-a_srs",proj,"-te",str(mygeo[0]),str(mygeo[1]),str(mygeo[2]),str(mygeo[3]),name1,name2])

def read_tree(treef):
    " Returns a pandas dataframe "

    A = pd.read_csv(treef,delimiter=' ',names=['linkno','start','end','dslink','uslink'])
    return A

def read_coord(coordf):
    " Returns a pandas dataframe "

    # we will only need columns contining lat and lon, there is an
    # empty column if twe read the entire file
    A = pd.read_csv(coordf,delimiter='\t',names=['lon','lat'],usecols=[1,2])
    return A

def haversine(point1, point2, miles=False):

    """
    Calculate the great-circle distance bewteen two points on the Earth surface.
    Uses Numpy functions

    """
    AVG_EARTH_RADIUS = 6371  # in km

    lat1, lng1 = point1
    lat2, lng2 = point2

    # convert all latitudes/longitudes from decimal degrees to radians
    lat1, lng1, lat2, lng2 = map(np.radians, (lat1, lng1, lat2, lng2))

    lat = lat2 - lat1
    lng = lng2 - lng1
    d = np.sin(lat*0.5)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(lng*0.5)**2
    h = 2*AVG_EARTH_RADIUS*np.arcsin(np.sqrt(d))
    if miles:
        return h * 0.621371  # in miles
    else:
        return h  # in kilometers

if __name__ == '__main__':
    getinflows(sys.argv[1:])
