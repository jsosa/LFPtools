#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# date: 23/dec/2017
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

def near_geo(ddsx,ddsy,XA):
    
    """
    Find nearest point in ddsx and ddsy np.array type arrays to a XA point
    XA should be defined as [x,y]. Nearest distance is calculated via Haversine
    formula (slower)
    """
    
    d = {'x':ddsx, 'y':ddsy, 'myx':XA[0], 'myy':XA[1]}
    df = pd.DataFrame(d)
    df['dis'] = df.apply(_pd_haversine,axis=1)
    ind = df['dis'].idxmin()
    dis = df.iloc[ind]['dis']
    return dis,ind

def near_euc(ddsx,ddsy,XA):
    
    """
    Find nearest point in ddsx and ddsy np.array type arrays to a XA point
    XA should be defined as [x,y]. Nearest distance is calculated via Euclidean
    formula (faster)
    """
    
    XA  = np.array([[XA[1],XA[0]]])
    XB  = np.vstack((ddsy,ddsx)).T
    dis = cdist(XA, XB, metric='euclidean').min()
    ind = cdist(XA, XB, metric='euclidean').argmin()
    return dis,ind

def neararray_geo(array,ddsx,ddsy,XA,tol):
    
    """
    Given an 2D array find nerest point to XA defined as [x,y]
    ddsx and ddsy are np.arrays containing information about row and
    columns. Mininum distance is calucalted using the Haversine formula
    """
    
    X,Y         = np.meshgrid(ddsx,ddsy)
    Xrav        = np.ravel(X)
    Yrav        = np.ravel(Y)
    dis,ind     = near_geo(Xrav,Yrav,XA)
    if dis<= tol:
        arr_y,arr_x = np.unravel_index(ind,array.shape)
        val         = array[arr_y,arr_x]
        x           = X[arr_y,arr_x]
        y           = Y[arr_y,arr_x]
        return x,y,val,dis
    else:
        return None        

def nearmask_geo(array,ddsx,ddsy,XA,tol):

    """
    Given an 2D array find nerest point to XA defined as [x,y]
    ddsx and ddsy are np.arrays containing information about row and
    columns. Mininum distance is calucalted using the Haversine formula
    It masks the array before
    """

    iy,ix       = np.where(array>0) # here is the mask!
    Xrav        = ddsx[ix]
    Yrav        = ddsy[iy]
    dis,ind     = near_geo(Xrav,Yrav,XA)
    if dis<= tol:
        val         = array[iy[ind],ix[ind]]
        x           = Xrav[ind]
        y           = Yrav[ind]
        return x,y,val,dis
    else:
        return None

def neararray_euc(array,ddsx,ddsy,XA,tol):
    
    """
    Given an 2D array find nerest point to XA defined as [x,y]
    ddsx and ddsy are np.arrays containing information about row and
    columns. Mininum distance is calucalted using the Euclidean formula
    """
    
    X,Y         = np.meshgrid(ddsx,ddsy)
    Xrav        = np.ravel(X)
    Yrav        = np.ravel(Y)
    dis,ind     = near_euc(Xrav,Yrav,XA)
    
    if dis<=tol:
        arr_y,arr_x = np.unravel_index(ind,array.shape)
        val         = array[arr_y,arr_x]
        x           = X[arr_y,arr_x]
        y           = Y[arr_y,arr_x]
        return x,y,val,dis
    else:
        return None        

def nearmask_euc(array,ddsx,ddsy,XA,tol):

    """
    Given an 2D array find nerest point to XA defined as [x,y]
    ddsx and ddsy are np.arrays containing information about row and
    columns. Mininum distance is calucalted using the Euclidean formula.
    It masks the array before
    """

    iy,ix       = np.where(array>0) # here is the mask!
    Xrav        = ddsx[ix]
    Yrav        = ddsy[iy]
    dis,ind     = near_euc(Xrav,Yrav,XA)
    if dis<= tol:
        val         = array[iy[ind],ix[ind]]
        x           = Xrav[ind]
        y           = Yrav[ind]
        return x,y,val,dis
    else:
        return None

def haversine(point1, point2, miles=False):

    """
    Calculate the great-circle distance bewteen two points on the Earth surface.
    Uses Numpy functions
    """
    
    AVG_EARTH_RADIUS = 6371  # in km

    lat1, lng1 = point1
    lat2, lng2 = point2

    # Convert all latitudes/longitudes from decimal degrees to radians
    lat1, lng1, lat2, lng2 = map(np.radians, (lat1, lng1, lat2, lng2))

    lat = lat2 - lat1
    lng = lng2 - lng1
    d = np.sin(lat*0.5)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(lng*0.5)**2
    h = 2*AVG_EARTH_RADIUS*np.arcsin(np.sqrt(d))
    if miles:
        return h * 0.621371  # in miles
    else:
        return h  # in kilometers

def _pd_haversine(row):

    """
    A wrapper to use Haversine formula in Pandas
    """
    
    point1 = [row['x'],row['y']]
    point2 = [row['myx'],row['myy']]
    return haversine(point1,point2)

def read_tree(treef):
    A = pd.read_csv(treef)
    A.set_index('index',inplace=True)
    return A

def read_coord(coordf):
    A = pd.read_csv(coordf)
    A.set_index('index',inplace=True)
    return A

def read_tree_taudem(treef):
    df = pd.read_csv(treef, sep='\t', names=['0','link_no','start_pnt','end_pnt','frst_ds'
                                             ,'frst_us','scnd_us','strahler','mon_pnt','shreve'])
    df.drop(['0'], axis=1, inplace=True)
    return df

def read_coord_taudem(coordf):
    df = pd.read_csv(coordf, sep='\t', names=['0','lon','lat','distance','elev','contr_area'])
    df.drop(['0'], axis=1, inplace=True)
    return df

def get_catchmentdir(filename):
    return os.path.dirname(filename) + '/'

def get_catchmentid(filename):
    return os.path.basename(os.path.dirname(filename))
