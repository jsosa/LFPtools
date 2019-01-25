#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

import os
import sys
import pdb
import getopt
import subprocess
import configparser
import numpy as np
import pandas as pd
import geopandas as gpd
import gdalutils as gu
from pyproj import Proj
from pyproj import transform
from shapely.geometry import Point


def getinflows_shell(argv):

    myhelp = '''
LFPtools v0.1

Name
----
getinflows

Description
-----------
Locate inflow points that can be used as boundary condition from a source e.g. discharge

Usage
-----
>> lfp-getinflows -i config.txt

Content in config.txt
---------------------
[getinflows]
ncf    = GeoTIFF file containing river network mask with means for example
ncproj = Projection mask e.g. epsg:3035
recf   = `Rec` file path
proj   = Projection output file e.g. epsg:4326
output = Output file
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

    ncf = str(config.get('getinflows', 'ncf'))
    ncproj = str(config.get('getinflows', 'ncproj'))
    recf = str(config.get('getinflows', 'recf'))
    proj = str(config.get('getinflows', 'proj'))
    output = str(config.get('getinflows', 'output'))

    getinflows(ncf, ncproj, recf, proj, output)


def getinflows(ncf, ncproj, recf, proj, output):

    print("    running getinflows.py...")

    fname = output

    # Reading XXX_rec.csv file
    rec = pd.read_csv(recf)

    # Creating inflow dataframe
    df_inf = pd.DataFrame()

    # Group by LINK
    recgrp = rec.groupby('link')
    for name, group in recgrp:

        mymean = []
        mydis = []
        df = group.copy()

        # Get lons and lats for every group
        lats = df['lat'].values
        lons = df['lon'].values

        for j in range(len(lons)):
            nearx1, neary1, ncmean, ncdis = find_nearest_mean_mask(
                ncf, ncproj, lons[j], lats[j], proj)
            if ncmean != None:
                mymean.append(ncmean)
                mydis.append(ncdis)
            else:
                mymean.append(np.nan)
                mydis.append(np.nan)
        df['mean'] = mymean
        df['dis'] = mydis
        df.dropna(inplace=True)

        # Finding best located points (close to JRC cell centers)
        df = df.loc[df.groupby('mean')['dis'].idxmin()]

        # Ensure points are sorted from larger to closest distance to the outlet
        df.sort_values(by='distance', ascending=False, inplace=True)

        # There will be links where no inflows are available
        # Then those links will be excluded and not considered
        if len(df) > 0:

            # Removing false inflow discharges
            flag = check_next_greater(df['mean'].values, 3)
            df['flag'] = flag
            df = df[df['flag'] == 1]

            # Removing last point in group, if there is one
            try:
                df.drop(index=df.index[0], inplace=True)
                df.drop(index=df.index[-1], inplace=True)
            except IndexError:
                pass

            # Updating inflows dataframe
            df_inf = pd.concat([df_inf, df])

    # Dropping duplicates from complete inflow list
    # df_inf.drop_duplicates('mean','first',inplace=True)

    # Dropping links with <2 points
    grouped = df_inf.groupby('link')
    df_new = pd.DataFrame()
    for name, group in grouped:
        if group['link'].size >= 2:
            df_new = pd.concat((df_new, group))

    # Final dataframe
    df_new = df_new.reset_index(drop=True)
    df_new.rename(columns={'lon': 'x', 'lat': 'y'}, inplace=True)

    # Create geodataframe
    gdf = gpd.GeoDataFrame(df_new, crs={'init': proj}, geometry=[
                           Point(xy) for xy in zip(df_new.x, df_new.y)])

    # Write geodataframe
    try:
        gdf.to_file(output, driver='GeoJSON')
    except:
        os.remove(output)
        gdf.to_file(output, driver='GeoJSON')


def find_nearest_mean_mask(ncf, ncproj, lon, lat, proj, thresh_mean=5, thresh_dis=2.5):
    """
    Apply a threshold to the mean discharge
    Based on the thresholded map, find nearest value to lon, lat in a given perimiter and var threhsold
    For JRC data default values are 5 m3s-1 and 2.5 Km
    """

    # JRC data set projection is EPSG:3035
    # It's required to convert to WGS84 to perform distance calculation
    crs_wgs84 = Proj(init=proj)
    crs_nc = Proj(init=ncproj)

    # Reading mean mask
    dat = gu.get_data(ncf)
    geo = gu.get_geo(ncf)

    # Create df, a pandas dataframe with values larger than "thresh_mean=5"
    df = gu.array_to_pandas(dat, geo, thresh_mean, 'ge')

    # Creating two new columns with projected values
    coords = transform(crs_nc, crs_wgs84, df['x'].values, df['y'].values)
    df['lon'] = coords[0]
    df['lat'] = coords[1]

    # Calcualte distance to lat and lon point to every point in the dataframe
    vec = gu.haversine.haversine_array(np.array(df['lat'], dtype='float32'),
                                       np.array(
        df['lon'], dtype='float32'),
        np.float32(lat),
        np.float32(lon))
    idx = np.argmin(vec)
    dis = gu.haversine.haversine(
        df.loc[idx, 'lat'], df.loc[idx, 'lon'], lat, lon)

    if dis <= thresh_dis:
        near_x = df.loc[idx, 'x']
        near_y = df.loc[idx, 'y']
        mymean = df.loc[idx, 'z']
        df = None
        return near_x, near_y, mymean, dis
    else:
        df = None
        return None, None, None, None


def check_next_greater(arr, thresh):
    """ Check for next greater point in an numpy array """

    flag = np.zeros(arr.shape[0], dtype=np.int)
    flag[0] = 1
    base_ = arr[0]
    for i in range(flag.size):
        next_ = arr[0+i]
        if (next_ > base_) & (next_/base_ <= thresh):
            base_ = next_
            flag[i] = 1
        else:
            j = 0
            while True:
                try:
                    next_2 = arr[i+j]
                    if (next_2 > base_) & (next_2/base_ <= thresh):
                        base_ = next_2
                        flag[i+j] = 1
                        break
                    else:
                        j += 1
                except IndexError:
                    break
    return flag


if __name__ == '__main__':
    getinflows_shell(sys.argv[1:])
