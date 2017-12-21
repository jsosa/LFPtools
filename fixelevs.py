#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# date: 21/apr/2017
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

import os
import sys
import getopt
import subprocess
import ConfigParser
import numpy as np
import shapefile
from osgeo import osr
import matplotlib.pyplot as plt
from itertools import groupby
from gdal_utils import *
from scipy.spatial.distance import cdist
from collections import defaultdict

import pdb # pdb.set_trace()

def fixelevs(argv):

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

    source = str(config.get('fixelevs','source'))
    output = str(config.get('fixelevs','output'))
    netf   = str(config.get('fixelevs','netf'))
    coordf = str(config.get('fixelevs','coordf'))
    treef  = str(config.get('fixelevs','treef'))
    proj   = str(config.get('fixelevs','proj'))
    method = str(config.get('fixelevs','method'))

    # the coord file correspond to the europe rivers, not the basin coord file!
    coord = np.genfromtxt(coordf, delimiter="\t")
    coord = np.delete(coord,0,1) # remove first column, is an empty column

    # the tree file correspond to the basin file
    tree = np.genfromtxt(treef, delimiter=" ")
    nlinks = len(open(treef,"r").read().split('\n')) - 1 # xxx_tre.txt file is generated with one a final empty ending line

    # database to fix
    elev   = np.array(shapefile.Reader(source).records(),dtype='float64')

    fname1 = output+"tmp"
    fname2 = output

    w = shapefile.Writer(shapefile.POINT)
    w.field('x')
    w.field('y')
    w.field('elevadj')

    # coordinates for bank elevations are based in river network mask
    net   = get_gdal_data(netf)
    geo   = get_gdal_geo(netf)
    iy,ix = np.where(net>0)
    xnet  = geo[8][ix]
    ynet  = geo[9][iy]

    for i in np.flipud(range(nlinks)):

        print "fixelevs.py - " + str(i)

        dem1 = []
        xx1  = []
        yy1  = []
        ii   = i

        while True:

            if nlinks == 1:
                link_beg = int(tree[1])
                link_end = int(tree[2]) + 1 # inclusive slicing
                linkds   = int(tree[3])
            else:
                link_beg = int(tree[ii,1])
                link_end = int(tree[ii,2]) + 1 # inclusive slicing
                linkds   = int(tree[ii,3])

            x1  = coord[link_beg:link_end,0]
            y1  = coord[link_beg:link_end,1]
            xx1 = np.append(xx1,x1)
            yy1 = np.append(yy1,y1)

            for ix1 in range(len(x1)):

                ielev = near(elev[:,1],elev[:,0],np.array([[x1[ix1],y1[ix1]]]))
                dem1  = np.append(dem1,elev[ielev,2])
            
            if nlinks  >  1: ii = np.where(tree[:,0]==linkds)
            if linkds == -1: break # the link is a downstream link?

        # calc bank elevation
        if method == 'yamazaki':
            adjusted_dem = bank4flood(dem1)
    
        for count in range(len(xx1)):
            x = xx1[count]
            y = yy1[count]
            inet = near(ynet,xnet,np.array([[x,y]])) # find the nearest point based on the river network raster mask
            w.point(xnet[inet],ynet[inet])
            w.record(xnet[inet],ynet[inet],adjusted_dem[count])

    w.save("%s.shp" % fname1)

    # write .prj file
    prj = open("%s.prj" % fname2, "w")
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj)
    prj.write(srs.ExportToWkt())
    prj.close()

    ### remove duplicated coordinates based on minium elevation value ###
    w = shapefile.Writer(shapefile.POINT)
    w.field('x')
    w.field('y')
    w.field('elevadj')
    elev  = np.array(shapefile.Reader(fname1+".shp").records(),dtype='float64')
    newx,newy,newz = check_duplicate_minimum(elev[:,0],elev[:,1],elev[:,2])
    for i in range(newx.size):
        if newx[i]!=-9999 and newy[i]!=-9999:
            w.point(newx[i],newy[i])
            w.record(newx[i],newy[i],newz[i])
    w.save("%s.shp" % fname2)
    
    nodata = -9999
    fmt    = "GTiff"
    name1  = output+".shp"
    name2  = output+".tif"
    subprocess.call(["gdal_rasterize","-a_nodata",str(nodata),"-of",fmt,"-tr",str(geo[6]),str(geo[7]),"-a","elevadj","-a_srs",proj,"-te",str(geo[0]),str(geo[1]),str(geo[2]),str(geo[3]),name1,name2])

def check_duplicate_minimum(x,y,z):

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
        ind  = np.argmin(z[D.values()[i]])
        indf = D.values()[i][ind]
        x[D.values()[i][ind]]=x1[D.values()[i][ind]]
        y[D.values()[i][ind]]=y1[D.values()[i][ind]]

    return x,y,z

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

    adjusted_dem=np.array(dem)

    for I in range(dem.size-1): # -1 to avoid error at boundary

        # bug on first and second elevation values
        if adjusted_dem[1]>adjusted_dem[0]:
            adjusted_dem[0]=adjusted_dem[1]

        if adjusted_dem[I+1]>adjusted_dem[I]:

            midind=I
            middem=adjusted_dem[midind]
            vecdem=adjusted_dem[midind:]

            # identify up pixel
            ii=0
            while vecdem[ii+1]>middem: # look downstream from pixel i and stop when pixel i+1 < midindex
                ii=ii+1
                if ii==vecdem.size-1: break # avoid problems at boundary downstream

            lastind=midind+ii+1

            zforw=adjusted_dem[midind:lastind]
            zsort=np.sort(zforw)
            
            zind=[] # indexes
            zmod=[] # adjusted elevation
            lmod=[] # cost function

            for J in range(zsort.size):

                # identify backward pixel
                jj=1
                while adjusted_dem[midind-jj]<=zsort[J]: # look backward from midindex and stop when pixel i-1 > mid-index
                    jj=jj+1
                    if jj>midind: break # avoid problems at boundary upstream

                backind=midind-jj+1

                z=adjusted_dem[backind:lastind] # extract DEM following backwardindex:forwardindex

                zind.append(range(backind,lastind))
                zmod.append(np.tile(zsort[J],(1,z.size))) # calc adjusted dem for every case
                lmod.append(np.sum(np.abs(z-zmod[J]))) # calc cost function for every case
            
            lmin=np.min(lmod)
            imin=np.where(lmod==lmin)[0][0]

            # print ""
            # print "    problem
            # print "    cost ....." + str(lmod[imin])
            # print "    indexes .." + str(np.float64(zind[imin]))
            # print "    before ..." + str(adjusted_dem[zind[imin]])
            # print "    after ...." + str(zmod[imin][0])
            # print ""

            adjusted_dem[zind[imin]]=zmod[imin] # final adjusted dem with minimum cost

    # remove flat banks
    adjusted_dem2 = avoid_flat_banks(adjusted_dem)

    # # DEBUG
    # fig, ax = plt.subplots()
    # ax.plot(range(dem.size),dem)
    # ax.plot(range(dem.size),adjusted_dem,'--')
    # plt.show()
    # # DEBUG

    return adjusted_dem2

def avoid_flat_banks(dem,thresh=0.001):

    for I in range(len(dem)-1):
        dem[I+1] = dem[I] - thresh

    return dem

def near(ddsx,ddsy,XA):

    XB  = np.vstack((ddsy,ddsx)).T
    dis = cdist(XA, XB, metric='euclidean').argmin()

    return dis

if __name__ == '__main__':
    fixelevs(sys.argv[1:])
