#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# date: 21/apr/2017
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

import sys
import getopt
import subprocess
import configparser
import numpy as np
import pandas as pd
import gdalutils
from lfptools import shapefile
from lfptools import misc_utils
from osgeo import osr
from sklearn import linear_model

def getslopes(argv):

    """
    Get slopes for the river network
    Finds nearestpoint in banks database to rec's coordinates
    """

    opts, args = getopt.getopt(argv,"i:")
    for o, a in opts:
        if o == "-i": inifile = a

    config = configparser.SafeConfigParser()
    config.read(inifile)

    source = str(config.get('getslopes','source'))
    output = str(config.get('getslopes','output'))
    netf   = str(config.get('getslopes','netf'))
    recf  = str(config.get('getslopes','recf'))
    proj   = str(config.get('getslopes','proj'))
    step   = int(config.get('getslopes','step'))

    print("    runnning get_slopes.py...")

    # Reading XXX_rec.csv file
    rec = pd.read_csv(recf)

    # Reading XXX_net.tif file
    geo = gdalutils.get_geo(netf)

    # Reading bank file (adjusted bank)
    elev = np.array(shapefile.Reader(source).records(),dtype='float64')

    # Initiate output shapefile
    w = shapefile.Writer(shapefile.POINT)
    w.field('x')
    w.field('y')
    w.field('slope')
    
    # Retrieving adjusted bank elevations from XXX_bnkfix.shp file
    # Values are stored in rec['bnk']
    bnkadj = []
    for i in rec.index:
        dis,ind = misc_utils.near_euc(elev[:,0],elev[:,1],(rec['lon'][i],
                                      rec['lat'][i]))
        bnkadj.append(elev[ind,2])
    rec['bnkadj'] = bnkadj

    # Calculating slopes
    # coordinates are grouped by REACH number
    rec['slopes'] = 0
    recgrp = rec.groupby('reach')
    for reach,df in recgrp:
        ids = df.index
        dem = df['bnkadj']
        # calc slopes
        slopes_vals = calc_slope_step(dem,df['lon'].values,df['lat'].values,step)
        rec['slopes'][ids] = slopes_vals

    # Writing .shp resulting file
    for i in rec.index:
        w.point(rec['lon'][i],rec['lat'][i])
        w.record(rec['lon'][i],rec['lat'][i],rec['slopes'][i])
    w.save("%s.shp" % output)

    # write .prj file
    prj = open("%s.prj" % output, "w")
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj)
    prj.write(srs.ExportToWkt())
    prj.close()
    
    # Writing .tif file
    nodata = -9999
    fmt    = "GTiff"
    name1  = output+".shp"
    name2  = output+".tif"
    subprocess.call(["gdal_rasterize","-a_nodata",str(nodata),"-of",fmt,"-tr",
                     str(geo[6]),str(geo[7]), "-a","slope","-a_srs",proj,"-te"
                     ,str(geo[0]),str(geo[1]),str(geo[2]),str(geo[3])
                     ,name1,name2])

def calc_slope_step(dem,x,y,step):

    myslp = np.ones(dem.size)*-9999
    
    # calculate distance by using haversine equation
    dis   = calc_dis_xy(x,y)

    # fit a linear regression by using Scikit learn, other more sophistcated
    # methods can be used to estimate the slope check on Linear regression methods
    # on Scikit learn website

    for i in range(len(dem)):

        left    = max(0,i-step)
        right   = min(len(dem),i+step+1) # +1 because inclusive slicing
        X_train = dis[left:right].reshape(-1,1)*1000 # *1000 -> to convert from kilometers in meters, reshape -> to requirement scikitlearn 
        Y_train = dem[left:right]
        regr    = linear_model.LinearRegression()
        regr.fit(X_train, Y_train)
        slp     = abs(regr.coef_)

        if slp <= 0.000001: slp = 0.0001
        myslp[i] = slp

        # # DEBUG DEBUG DEBUG
        # plt.scatter(X_train, Y_train,  color='black')
        # plt.scatter(X_train[step],Y_train[step], color='red')
        # plt.plot(X_train, regr.predict(X_train), color='blue', linewidth=3)
        # plt.show()

    return myslp

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

def calc_dis_xy(x,y):
    dis = np.zeros(x.size)
    for i in range(len(dis)):
        if i > 0:
            dis[i] = haversine([y[i],x[i]],[y[i-1],x[i-1]])
        discum = np.cumsum(dis)
    return discum

if __name__ == '__main__':
    getslopes(sys.argv[1:])
