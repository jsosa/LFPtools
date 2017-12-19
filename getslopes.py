#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# date: 21/apr/2017
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

import os,sys,getopt,subprocess
import ConfigParser
import numpy as np
import shapefile
import matplotlib.pyplot as plt
from osgeo import osr
from gdal_utils import *
from scipy.spatial.distance import cdist
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
from collections import defaultdict

import pdb
# pdb.set_trace()

def getslopes(argv):

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

    opts, args = getopt.getopt(argv,"i:")
    for o, a in opts:
        if o == "-i": inifile = a

    config = ConfigParser.SafeConfigParser()
    config.read(inifile)

    source = str(config.get('getslopes','source'))
    output = str(config.get('getslopes','output'))
    netf   = str(config.get('getslopes','netf'))
    coordf = str(config.get('getslopes','coordf'))
    treef  = str(config.get('getslopes','treef'))
    proj   = str(config.get('getslopes','proj'))
    method = str(config.get('getslopes','method'))
    
    try:
        # required if linear_regr_step is used!
        step   = int(config.get('getslopes','step'))
    except: pass

    # the coord file correspond to the europe rivers, not the basin coord file!
    coord = np.genfromtxt(coordf, delimiter="\t")
    coord = np.delete(coord,0,1) # remove first column, is an empty column

    # the tree file correspond to the basin file
    tree = np.genfromtxt(treef, delimiter=" ")
    nlinks = len(open(treef,"r").read().split('\n')) - 1 # xxx_tre.txt file is generated with one a final empty ending line

    # database to estimate slopes, it should be the bank .shp file
    elev   = np.array(shapefile.Reader(source).records(),dtype='float64')

    fname1 = output+"tmp"
    fname2 = output

    w = shapefile.Writer(shapefile.POINT)
    w.field('x')
    w.field('y')
    w.field('slope')

    # coordinates for bank elevations are based in river network mask
    net   = get_gdal_data(netf)
    geo   = get_gdal_geo(netf)
    iy,ix = np.where(net>0)
    xnet  = geo[8][ix]
    ynet  = geo[9][iy]

    for i in np.flipud(range(nlinks)):

        print "getslopes.py - " + str(i)

        dem1 = []

        if nlinks == 1:
            link_beg = int(tree[1])
            link_end = int(tree[2]) + 1 # inclusive slicing
        else:
            link_beg = int(tree[i,1])
            link_end = int(tree[i,2]) + 1 # inclusive slicing

        x1  = coord[link_beg:link_end,0]
        y1  = coord[link_beg:link_end,1]

        for ix1 in range(len(x1)):
            ielev = near(elev[:,1],elev[:,0],np.array([[x1[ix1],y1[ix1]]]))
            dem1  = np.append(dem1,elev[ielev,2])

        # calc slope
        if method == 'classic':
            myslope = calc_slope(dem1,x1,y1,1000)

        elif method == 'linear_regr_step':
            myslope = calc_slope_step(dem1,x1,y1,step)

        for count in range(len(myslope)):
            x = x1[count]
            y = y1[count]
            inet = near(ynet,xnet,np.array([[x,y]])) # find the nearest point based on the river network raster mask
            w.point(xnet[inet],ynet[inet])
            w.record(xnet[inet],ynet[inet],myslope[count])

    w.save("%s.shp" % fname1)

    # write .prj file
    prj = open("%s.prj" % fname1, "w")
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj)
    prj.write(srs.ExportToWkt())
    prj.close()

    ### remove duplicated coordinates based on maximum slope value ###
    w = shapefile.Writer(shapefile.POINT)
    w.field('x')
    w.field('y')
    w.field('slope')
    shp = np.array(shapefile.Reader(fname1+".shp").records(),dtype='float64')
    newx,newy,newz = check_duplicate_max(shp[:,0],shp[:,1],shp[:,2])
    for i in range(newx.size):
        if newx[i]!=-9999 and newy[i]!=-9999:
            w.point(newx[i],newy[i])
            w.record(newx[i],newy[i],newz[i])
    w.save("%s.shp" % fname2)
    
    nodata = -9999
    fmt    = "GTiff"
    name1  = output+".shp"
    name2  = output+".tif"
    subprocess.call(["gdal_rasterize","-a_nodata",str(nodata),"-of",fmt,"-tr",str(geo[6]),str(geo[7]),"-a","slope","-a_srs",proj,"-te",str(geo[0]),str(geo[1]),str(geo[2]),str(geo[3]),name1,name2])

def calc_slope(dem,x,y,res):

    myslp = np.ones(dem.size)*-9999

    for i in range(len(dem)):

        if i == len(dem)-1:
            dis = haversine([y[i],x[i]],[y[i-1],x[i-1]])
            slp = abs(dem[i]-dem[i-1])/(dis*res) # units are [m/m]
        else:
            dis = haversine([y[i],x[i]],[y[i+1],x[i+1]])
            slp = abs(dem[i]-dem[i+1])/(dis*res) # units are [m/m]

        if slp == 0: slp = 0.0001

        myslp[i] = slp

    return myslp

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
        X_train = dis[left:right].reshape(-1,1)*1000 # *1000 -> to convert to kilometers in meters, reshape -> to requirement scikitlearn 
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

def near(ddsx,ddsy,XA):

    XB  = np.vstack((ddsy,ddsx)).T
    ind = cdist(XA, XB, metric='euclidean').argmin()

    return ind

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

def check_duplicate_max(x,y,z):

    """
    After applying Yamzaki method several coordinates are duplicated
    elevation values varies because the method is applied at each reach
    then, is necesary to take the minimum value!, maybe?

    """
    x1 = np.copy(x)
    y1 = np.copy(y)
    mylist = zip(x,y)
    D = defaultdict(list)
    for i,item in enumerate(mylist):
        D[item].append(i)
    D = {k:v for k,v in D.items() if len(v)>1}

    for i in range(len(D)):
        x[D.values()[i]]=-9999
        y[D.values()[i]]=-9999
        
        # minimum is value is taken, can be changed!
        ind  = np.argmax(z[D.values()[i]])
        indf = D.values()[i][ind]
        x[D.values()[i][ind]]=x1[D.values()[i][ind]]
        y[D.values()[i][ind]]=y1[D.values()[i][ind]]

    return x,y,z

if __name__ == '__main__':
    getslopes(sys.argv[1:])
