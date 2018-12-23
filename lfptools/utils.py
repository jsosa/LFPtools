#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: j.sosa@bristol.ac.uk / sosa.jeison@gmail.com

import os
import shutil
from glob import glob
import gzip
import numpy as np
import pandas as pd
import gdalutils
from osgeo import osr


def _secs_to_time(df, date1):

    df = df.copy()
    conversion = 86400  # 86400s = 1day
    df['time'] = pd.to_datetime(
        df['Time']/conversion, unit='D', origin=pd.Timestamp(date1))
    df.set_index(df['time'], inplace=True)
    del df['Time']
    del df['time']
    return df


def _hours_to_time(df, date1):

    df = df.copy()
    conversion = 24  # 24h = 1day
    df['time'] = pd.to_datetime(
        df['Time']/conversion, unit='D', origin=pd.Timestamp(date1))
    df.set_index(df['time'], inplace=True)
    del df['Time']
    del df['time']
    return df


def _get_lineno(filename, phrase):
    with open(filename, 'r') as f:
        for num, line in enumerate(f):
            if phrase in line:
                return num


def read_mass(filename, date1='1990-01-01'):

    df = pd.read_csv(filename, delim_whitespace=True)
    df = _secs_to_time(df, date1)
    df['res'] = np.arange(0, df.index.size)
    return df


def read_discharge(filename, date1='1990-01-01'):

    line = _get_lineno(filename, 'Time') + 1  # inclusive slicing
    df = pd.read_csv(filename, skiprows=range(0, line),
                     header=None, delim_whitespace=True)
    df.rename(columns={0: 'Time'}, inplace=True)
    df = _secs_to_time(df, date1)
    return df


def read_stage(filename, date1='1990-01-01'):

    line = _get_lineno(filename, 'Time') + 1  # inclusive slicing
    df = pd.read_csv(filename, skiprows=range(0, line),
                     header=None, delim_whitespace=True)
    df.rename(columns={0: 'Time'}, inplace=True)
    df = _secs_to_time(df, date1)
    return df


def read_stage_locs(filename):

    str_line = _get_lineno(filename, 'Stage information') + 1
    end_line = _get_lineno(filename, 'Output, depths:') - 1
    df = pd.read_csv(filename, header=None, delim_whitespace=True,
                     skiprows=range(0, str_line), nrows=end_line-str_line,
                     index_col=0, names=['x', 'y', 'elev'])
    return df


def read_bci(filename):

    return pd.read_csv(filename, skiprows=1, delim_whitespace=True,
                       names=['boundary', 'x', 'y', 'type', 'name'])


def read_bdy(filename, bcifile, date1='1990-01-01'):

    phrase = 'hours'
    bdy = pd.DataFrame()
    with open(filename, 'r') as f:
        for num, line in enumerate(f):
            if phrase in line:
                start = num + 1
                lines = int(line.split(' ')[0])
                total = start + lines
                df = pd.read_csv(filename, skiprows=start, nrows=total-start,
                                 header=None, delim_whitespace=True)
                bdy = pd.concat([bdy, df[0]], axis=1)
    bdy['Time'] = df[1]
    bdy = _hours_to_time(bdy, date1)
    bdy.columns = read_bci(bcifile).name.values

    return bdy


def read_par(filename):

    df = pd.read_csv(filename, delim_whitespace=True, header=None, index_col=0)
    df.fillna('', inplace=True)
    df.columns = ['']
    df.index.name = ''
    return df


def get_ascii_geo(filename, proj4=None):
    """
    Reads LISFLOOD-FP outputs either in gz or wd file extension GEO info
    Assumes WGS_84 projection by default
    """

    ext = os.path.basename(filename).split('.')[-1]

    if ext == 'gz':
        file = _uncompress_gz(filename)
        geo = gdalutils.get_geo(file, proj4=proj4)
        os.remove(file)
    else:
        geo = gdalutils.get_geo(filename, proj4=proj4)

    return geo


def get_ascii_dat(filename):
    """
    Reads LISFLOOD-FP outputs either in gz or wd file extension DAT info
    """

    ext = os.path.basename(filename).split('.')[-1]

    if ext == 'gz':
        file = _uncompress_gz(filename)
        data = gdalutils.get_data(file)
        os.remove(file)
    else:
        data = gdalutils.get_data(filename)

    return data
