#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# date: 13/dec/2017
# mail: j.sosa@bristol.ac.uk / sosa.jeison@gmail.com

cimport cython
import numpy as np
cimport numpy as np
from libc.math cimport sin,asin,cos,pow,sqrt,M_PI

def cy_rastermask(np.float64_t[:,:] data, np.int16_t[:,:] mask):

    """
    Mask a raster based on a predefined mask
    """

    cdef np.int32_t M = data.shape[0]
    cdef np.int32_t N = data.shape[1]
    cdef np.int32_t n,m

    for m in range(M):
        for n in range(N):
            if mask[m,n] >= 1:
                data[m,n] = data[m,n]
            elif mask[m,n] == 0:
                data[m,n] = 0
    return data

def cy_rasterthreshold(np.float64_t[:,:] data, np.float64_t thresh, np.float64_t nodata):

    """
    Threshold a raster ex. accumulation raster to get river network mask (1,0)
    """

    cdef np.int32_t M = data.shape[0]
    cdef np.int32_t N = data.shape[1]
    cdef np.int32_t n,m

    for m in range(M):
        for n in range(N):
            if data[m,n] >= thresh:
                data[m,n] = 1
            elif data[m,n] == nodata:
                data[m,n] = data[m,n]
            else:
              data[m,n] = 0
    return data

def cy_directions_tau(np.int16_t[:,:] data, np.int16_t nodata):

    """
    Function to use in Shell to change convetion from a DIR file
    HydroSHEDS uses ESRI convention 128,64,32,.. this script
    changes these numbers to TauDEM convention, 1,2,3...
    """

    cdef np.int32_t M = data.shape[0]
    cdef np.int32_t N = data.shape[1]
    cdef np.int32_t n,m

    for m in range(M):
        for n in range(N):
            if data[m,n] == 8:
                data[m,n] = 6
            elif data[m,n] == 2:
                data[m,n] = 8
            elif data[m,n] == 128:
                data[m,n] = 2
            elif data[m,n] == 4:
                data[m,n] = 7
            elif data[m,n] == 32:
                data[m,n] = 4
            elif data[m,n] == 16:
                data[m,n] = 5
            elif data[m,n] == 64:
                data[m,n] = 3
            elif data[m,n] == 0:
                data[m,n] = nodata
            elif data[m,n] == 247:
                data[m,n] = nodata
            elif data[m,n] == 255:
                data[m,n] = nodata
    return data

def cy_d82d4(np.int16_t[:,:] data, np.int16_t nodata):

    """
    Returns direction and river netwrok maps in D4
    """

    cdef np.int32_t M = data.shape[0]
    cdef np.int32_t N = data.shape[1]
    cdef np.int32_t n,m

    for m in range(M):
        for n in range(N):
          if data[m,n] == 8:
            data[m,n] = 7
            try:
              data[m+1,n] = 1
            except:
              pass

          elif data[m,n] == 6:
            data[m,n] = 7
            try:
              data[m+1,n] = 5
            except:
              pass

          elif data[m,n] == 4:
            data[m,n] = 3
            try:
              data[m-1,n] = 5
            except:
              pass

          elif data[m,n] == 2:
            data[m,n] = 3
            try:
              data[m-1,n] = 1
            except:
              pass

    cdef np.int16_t[:,:] net = np.ones((M,N),dtype=np.int16) * np.greater(data,0)

    return (data,net)


@cython.wraparound(False)
@cython.boundscheck(False)
def calc_area(np.int16_t nx, np.int16_t ny, np.float32_t resx, np.float32_t resy, np.float32_t[:] x, np.float32_t[:] y):

    cdef np.int16_t i,j
    cdef np.float32_t xx,yy,x1,x2,y1,y2
    cdef np.float32_t[:,:] area = np.zeros((ny,nx),dtype=np.float32)

    for j in range(ny):
        for i in range(nx):
            x1 = x[i] - resx
            x2 = x[i] + resx
            y1 = y[j] + resy
            y2 = y[j] - resy
            xx = haversine(y1,x1,y1,x2)
            yy = haversine(y1,x1,y2,x1)
            area[j,i] = xx*yy

    return area

@cython.wraparound(False)
@cython.boundscheck(False)
def haversine(np.float32_t lat1, np.float32_t lng1, np.float32_t lat2, np.float32_t lng2):
    
    cdef np.int32_t AVG_EARTH_RADIUS = 6371  # in km
    cdef np.float32_t lat,lng,d,h

    # Convert all latitudes/longitudes from decimal degrees to radians
    lat1 = lat1*M_PI/180
    lng1 = lng1*M_PI/180
    lat2 = lat2*M_PI/180
    lng2 = lng2*M_PI/180

    lat = lat2 - lat1
    lng = lng2 - lng1
    d = pow(sin(lat*0.5),2) + cos(lat1)*cos(lat2)*pow(sin(lng*0.5),2)
    h = 2*AVG_EARTH_RADIUS*asin(sqrt(d))
    
    return h
