#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: j.sosa@bristol.ac.uk / sosa.jeison@gmail.com

import os
import shutil
import gzip
import numpy as np
import pandas as pd
import gdalutils
from osgeo import osr


# TODO:
    # update listing dir by using from glob import glob
    # include projection info when creating the netcdf file


def _secs_to_time(df,date1):

    conversion = 86400 # 86400s = 1day
    df['time'] = pd.to_datetime(df['Time']/conversion,unit='D',origin=pd.Timestamp(date1)) 
    df.set_index(df['time'],inplace=True)
    del df['Time']
    del df['time']
    return df


def _hours_to_time(df,date1):

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


def _stack_variable(path,ext):
    mylis = []
    for file in os.listdir(path):
        if file.endswith(ext):
            fname = os.path.join(path,file)
            mydat = get_ascii_data(fname)
            mylis.append(mydat)
    myvar = np.dstack(mylis)

    # Replace zeros by np.nan
    myvar[myvar==0] = np.nan

    # Retrieve geo information from last file
    mygeo = get_ascii_geo(fname)

    return myvar,mygeo


def results_to_nc(path,mass,date1='1990-01-01'):

    import xarray as xr

    # Mass file is used only to read output dates
    t = pd.read_csv(mass, delim_whitespace=True, usecols=[0,0])
    t = _secs_to_time(t,date1)

    wd_dat,wd_geo = _stack_variable(path,'.wd.gz')
    wdfp_dat,wdfp_geo = _stack_variable(path,'.wdfp.gz')
    elev_dat,elev_geo = _stack_variable(path,'.elev.gz')
    Fwidth_dat,Fwidth_geo = _stack_variable(path,'.Fwidth.gz')
    Qx_dat,Qx_geo = _stack_variable(path,'.Qx.gz')
    Qy_dat,Qy_geo = _stack_variable(path,'.Qy.gz')
    Qcx_dat,Qcx_geo = _stack_variable(path,'.Qcx.gz')
    Qcy_dat,Qcy_geo = _stack_variable(path,'.Qcy.gz')

    # Creating xarray DataArray to be exported in NetCDF
    x     = wd_geo[8]
    y     = wd_geo[9]
    time  = t.index.values
    
    foo                            = xr.Dataset()
    foo['wd']                      = (('y','x','time'), wd_dat)
    foo['wdfp']                    = (('y','x','time'), wdfp_dat)
    foo['elev']                    = (('y','x','time'), elev_dat)
    foo['Fwidth']                  = (('y','x','time'), Fwidth_dat)
    # foo['Qx']                      = (('y','x','time'), Qx_dat)
    # foo['Qy']                      = (('y','x','time'), Qy_dat)
    # foo['Qcx']                     = (('y','x','time'), Qcx_dat)
    # foo['Qcy']                     = (('y','x','time'), Qcy_dat)
    # foo['Q']                       = xr.ufuncs.sqrt(foo.Qx**2 + foo.Qy**2)
    # foo['Qc']                      = xr.ufuncs.sqrt(foo.Qy**2 + foo.Qcy**2)
    foo.coords['x']                = (('x'), x)
    foo.coords['y']                = (('y'), y)
    foo.coords['time']             = time
    foo.attrs['creator_name']      = 'Jeison Sosa'
    foo.attrs['creator_institute'] = 'University of Bristol'
    foo.attrs['creator_email']     = 'j.sosa@bristol.ac.uk'

    foo.to_netcdf(path+"res.nc",encoding={'wd': {'zlib': True}, 'wdfp': {'zlib': True}})
 

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
    COnverts ASCII file into GTiff
    """

    myarray = gdalutils.get_data(filename)
    mygeo = gdalutils.get_geo(filename,proj4=proj4)

    return gdalutils.write_raster(myarray,output,mygeo,fmt,mygeo[11])
