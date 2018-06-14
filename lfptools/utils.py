#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: j.sosa@bristol.ac.uk / sosa.jeison@gmail.com

import os
import shutil
from glob import glob
import gzip
import xarray as xr
import numpy as np
import pandas as pd
import gdalutils
from osgeo import osr
import gdalutils.extras.haversine as haversine


# TODO:
    # include projection info when creating the netcdf file


def _secs_to_time(df,date1):

    df = df.copy()
    conversion = 86400 # 86400s = 1day
    df['time'] = pd.to_datetime(df['Time']/conversion,unit='D',origin=pd.Timestamp(date1)) 
    df.set_index(df['time'],inplace=True)
    del df['Time']
    del df['time']
    return df


def _hours_to_time(df,date1):

    df = df.copy()
    conversion = 24 # 24h = 1day
    df['time'] = pd.to_datetime(df['Time']/conversion,unit='D',origin=pd.Timestamp(date1)) 
    df.set_index(df['time'],inplace=True)
    del df['Time']
    del df['time']
    return df


def _get_lineno(filename,phrase):
    with open(filename,'r') as f:
        for num, line in enumerate(f):
            if phrase in line:
                return num
            

def read_mass(filename,date1='1990-01-01'):

    df = pd.read_csv(filename, delim_whitespace=True)
    df = _secs_to_time(df,date1)
    df['res'] = np.arange(0,df.index.size)
    return df


def read_discharge(filename,date1='1990-01-01'):

    line = _get_lineno(filename,'Time') + 1 # inclusive slicing
    df = pd.read_csv(filename,skiprows=range(0,line),header=None,delim_whitespace=True)
    df.rename(columns={0:'Time'},inplace=True)
    df = _secs_to_time(df,date1)
    return df


def read_stage(filename,date1='1990-01-01'):

    line = _get_lineno(filename,'Time') + 1 # inclusive slicing
    df = pd.read_csv(filename,skiprows=range(0,line),header=None,delim_whitespace=True)
    df.rename(columns={0:'Time'},inplace=True)
    df = _secs_to_time(df,date1)
    return df


def read_stage_locs(filename):
    
    str_line = _get_lineno(filename,'Stage information') +1
    end_line = _get_lineno(filename,'Output, depths:') - 1
    df = pd.read_csv(filename,header=None,delim_whitespace=True,
                     skiprows=range(0,str_line),nrows=end_line-str_line,
                     index_col=0,names=['x','y','elev'])
    return df


def read_bci(filename):

    return pd.read_csv(filename,skiprows=1,delim_whitespace=True,
                     names=['boundary','x','y','type','name'])


def read_bdy(filename,bcifile,date1='1990-01-01'):

    phrase = 'hours'
    bdy = pd.DataFrame()
    with open(filename,'r') as f:
        for num, line in enumerate(f):
            if phrase in line:
                start = num + 1
                lines = int(line.split(' ')[0])
                total = start+ lines
                df = pd.read_csv(filename,skiprows=start,nrows=total-start,
                                 header=None,delim_whitespace=True)
                bdy = pd.concat([bdy,df[0]],axis=1)
    bdy['Time'] = df[1]
    bdy = _hours_to_time(bdy,date1)
    bdy.columns = read_bci(bcifile).name.values

    return bdy


def read_par(filename):

    df = pd.read_csv(filename, delim_whitespace=True, header=None, index_col=0)
    df.fillna('',inplace=True)
    df.columns=['']
    df.index.name=''
    return df


def _stack_variable(path,ext,compress=''):

    files = glob(path+'*.'+ext+compress)
    lengt = len(files)
    mylis = []
    for i in range(lengt):
        num   = '%04d' % i 
        fname = path + 'res' + '-' + num + '.' + ext+compress
        mydat = get_ascii_dat(fname)
        mylis.append(mydat)

    myvar = np.dstack(mylis)

    # Replace zeros by np.nan
    myvar[myvar==0] = np.nan

    # Retrieve geo information from last file
    mygeo = get_ascii_geo(fname)

    return myvar,mygeo


def results_to_nc(path,mass,var,compress,outf,date1='1990-01-01'):

    if compress == True:
        compress = '.gz'
    else:
        compress = ''

    # Mass file is used only to read output dates
    t = pd.read_csv(mass, delim_whitespace=True, usecols=[0,0])
    t = _secs_to_time(t,date1)

    dat,geo = _stack_variable(path,var,compress)

    x    = geo[8]
    y    = geo[9]
    time = t.index.values
    
    foo                            = xr.Dataset()
    foo[var]                       = (('y','x','time'), dat)
    foo.coords['x']                = (('x'), x)
    foo.coords['y']                = (('y'), y)
    foo.coords['time']             = time
    foo.attrs['creator_name']      = 'Jeison Sosa'
    foo.attrs['creator_institute'] = 'University of Bristol'
    foo.attrs['creator_email']     = 'j.sosa@bristol.ac.uk'

    foo.to_netcdf(outf,encoding={var: {'zlib': True}})
 

def _uncompress_gz(filename):

    path = os.path.dirname(filename) + '/'
    name = 'temp.txt'
    file = path + name

    with gzip.open(filename, 'rb') as f_in:
        with open(file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    return file


def get_ascii_geo(filename,proj4=None):

    """
    Reads LISFLOOD-FP outputs either in gz or wd file extension GEO info
    Assumes WGS_84 projection by default
    """

    ext = os.path.basename(filename).split('.')[-1]

    if ext == 'gz':
        file = _uncompress_gz(filename)
        geo  = gdalutils.get_geo(file,proj4=proj4)
        os.remove(file)
    else:
        geo = gdalutils.get_geo(filename,proj4=proj4)

    return geo


def get_ascii_dat(filename):

    """
    Reads LISFLOOD-FP outputs either in gz or wd file extension DAT info
    """
    
    ext = os.path.basename(filename).split('.')[-1]
    
    if ext == 'gz':
        file = _uncompress_gz(filename)
        data  = gdalutils.get_data(file)
        os.remove(file)
    else:
        data = gdalutils.get_data(filename)

    return data


def to_tif(filename,output,proj4=None,fmt="Float64"):

    """
    Converts ASCII file into GTiff
    """

    myarray = gdalutils.get_data(filename)
    mygeo   = gdalutils.get_geo(filename,proj4=proj4)

    return gdalutils.write_raster(myarray,output,mygeo,fmt,mygeo[11])
