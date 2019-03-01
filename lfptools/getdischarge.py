#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: j.sosa@bristol.ac.uk / sosa.jeison@gmail.com

import sys
import getopt
import configparser
import numpy as np
import pandas as pd
import xarray as xr
import gdalutils as gu
import geopandas as gpd
from pyproj import Proj
from pyproj import transform


def getdischarge_shell(argv):

    myhelp = '''
LFPtools v0.1

Name
----
getdischarge

Description
-----------
Retrieve discharge from a netCDF source following a predefined inflows locations (e.g. from lfp-getinflows)

Usage
-----
>> lfp-getdischarge -i config.txt

Content in config.txt
---------------------
[getdischarge]
ncf    = NetCDF source file
ncproj = Projection netCDF file e.g. epsg:3035
infshp = Inflow locations e.g. from lfp-getinflows
proj   = Output projection e.g. epsg:4326
output = Output file, a CSV
date1  = Start date e.g. 1990-01-01
date2  = End date e.g. 1990-01-30
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

    ncf = str(config.get('getdischarge', 'ncf'))
    ncproj = str(config.get('getdischarge', 'ncproj'))
    ncxlabel = str(config.get('getdischarge', 'ncxlabel'))
    ncylabel = str(config.get('getdischarge', 'ncylabel'))
    ncdatlbl = str(config.get('getdischarge', 'ncdatlbl'))
    infshp = str(config.get('getdischarge', 'infshp'))
    proj = str(config.get('getdischarge', 'proj'))
    output = str(config.get('getdischarge', 'output'))
    date1 = str(config.get('getdischarge', 'date1'))
    date2 = str(config.get('getdischarge', 'date2'))

    getdischarge(ncf, ncproj, ncxlabel, ncylabel, ncdatlbl,
                 infshp, proj, output, date1, date2)


def getdischarge(ncf, ncproj, ncxlabel, ncylabel, ncdatlbl, infshp, proj, output, date1, date2):

    print("    running getdischarge.py...")

    # Reading inflows
    gdf = gpd.read_file(infshp)

    # Getting nearest discharge point
    def _pd_find_nearest(row):
        x = row.x
        y = row.y
        near_x, near_y = find_nearest(
            ncf, ncproj, ncxlabel, ncylabel, x, y, proj)
        row['near_x'] = near_x
        row['near_y'] = near_y
        return row

    df = gdf.apply(_pd_find_nearest, axis=1)

    # Retriving discharges for the time specified and specified inflow points
    res = pd.DataFrame()
    for i in range(len(df)):
        near_x = df.iloc[i]['near_x']
        near_y = df.iloc[i]['near_y']
        tserie = get_data(ncf, ncdatlbl, ncxlabel, ncylabel,
                          near_x, near_y, date1, date2)
        dis = tserie
        dis['index'] = i
        dis['columns'] = dis.index.astype(str)
        dispivot = dis.pivot(
            columns='columns', index='index', values='discharge')
        res = pd.concat([res, dispivot])

    df1 = pd.concat([df, res], axis=1)

    # Writing CSV file
    df1.to_csv(output)


def get_data(ncf, ncdatlbl, ncxlabel, ncylabel, x, y, date1='1990-01-01', date2='2014-12-31'):
    """
    Retrieve array based on nearest x,y coordinates

    x : longitude same projection as source
    y : latitude same projection as source
    """

    dat = xr.open_dataset(ncf)
    mytim = dat.sel(time=slice(date1, date2))
    mydis = mytim.sel({ncxlabel: x, ncylabel: y}, method="nearest")
    df = mydis[ncdatlbl].to_pandas().to_frame()
    df.columns = ['discharge']
    return df


def find_nearest(ncf, ncproj, ncxlabel, ncylabel, lon, lat, proj):
    """
    Find nearest point in discharge dataset based on inflow points

    ncf : netcdf file
    ncroj : projection e.g. epsg:3035
    lon : longitude
    lat : latitude
    proj : lon and lat projection e.g. epsg:4326
    """

    crs_wgs84 = Proj(init=proj)
    crs_nc = Proj(init=ncproj)

    # Reading netcdf file
    dat = xr.open_dataset(ncf)

    # Transforming between projections
    x, y = transform(crs_wgs84, crs_nc, lon, lat)

    # Retrieve near x and y
    near = dat.sel({ncxlabel: x, ncylabel: y}, method="nearest")
    near_x = np.float64(near[ncxlabel].values)
    near_y = np.float64(near[ncylabel].values)

    return near_x, near_y


if __name__ == '__main__':
    getdischarge_shell(sys.argv[1:])
