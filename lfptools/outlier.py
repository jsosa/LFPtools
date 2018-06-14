#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

import sys
import getopt
import subprocess
import configparser
import numpy as np
import pandas as pd
from lfptools.fixelevs import bank4flood
import gdalutils
import gdalutils.extras.shapefile as shapefile
import gdalutils.extras.haversine as haversine
from osgeo import osr

def main(argv):

    # This program was made to fix or apply consistency in datasets. For example,
    # at some points depth can be 30 m surrounded of values of 3 m and 5 m, clearly
    # this is an inconsistency in the reach. This value is removed and replaces by
    # an interpolated value get it using correct values.

    # It creates outputs in shapefile and geotiff formats.

    opts, args = getopt.getopt(argv,"i:")
    for o, a in opts:
        if o == "-i": inifile  = a

    config = configparser.SafeConfigParser()
    config.read(inifile)

    recf = str(config.get('outlier','recf'))
    proj = str(config.get('outlier','proj'))
    netf = str(config.get('outlier','netf'))
    shpf = str(config.get('outlier','shpf'))
    labl = str(config.get('outlier','labl'))
    outf = str(config.get('outlier','outf'))
    meth = str(config.get('outlier','meth'))

    print("    running outlier.py...")

    outlier(recf,proj,netf,shpf,labl,outf,meth)

def outlier(recf,proj,netf,shpf,labl,outf,method):

    # Loading xxx_rec file
    rec = pd.read_csv(recf)

    # Loading shapefile
    df  = gdalutils.shp_to_pandas(shpf)

    # Assign values in rec file
    rec = assign_val(rec,df,labl)

    # Remove outliers in new rec file
    if method == 'qua':
        rec = remove_outliers_qua(rec,'link',labl)
    elif method == 'mad':
        rec = remove_outliers_mad(rec,'link',labl)
    elif method == 'yama':
        rec = remove_outliers_yama(rec,'reach',labl)
    else:
        sys.exit('ERROR Method not recognized')

    # Saving shapefile and tif files
    save_shp_tif(netf,outf,proj,rec,labl)

def manualfix(recf,proj,netf,shpf,labl,outf,x,y,val):

    # Loading xxx_rec file
    rec = pd.read_csv(recf)

    # Loading shapefile
    df  = gdalutils.shp_to_pandas(shpf)

    # Assign values in rec file
    rec = assign_val(rec,df,labl)

    for i in range(len(x)):
        dis = haversine.haversine_array(np.array(rec['lat'],dtype='float32'),
                                        np.array(rec['lon'],dtype='float32'),
                                        np.float32(y[i]),
                                        np.float32(x[i]))
        idx = np.argmin(dis)
        rec.loc[idx,labl] = val[i]

    # Saving shapefile and tif files
    save_shp_tif(netf,outf,proj,rec,labl)

def save_shp_tif(netf,outf,proj,rec,labl):

    # Save rec in shapefile
    w = shapefile.Writer(shapefile.POINT)
    w.field('x')
    w.field('y')
    w.field(labl)

    # Writing .shp resulting file
    for x,y,value in zip(rec['lon'],rec['lat'],rec[labl]):
        w.point(x,y)
        w.record(x,y,value)
    w.save("%s.shp" % outf)

    # write .prj file
    prj = open("%s.prj" % outf, "w")
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj)
    prj.write(srs.ExportToWkt())
    prj.close()

    geo = gdalutils.get_geo(netf)

    fmt    = "GTiff"
    nodata = -9999
    name1  = outf+".shp"
    name2  = outf+".tif"
    subprocess.call(["gdal_rasterize","-a_nodata",str(nodata),
                                      "-of",fmt,
                                      "-tr",str(geo[6]),str(geo[7]),
                                      "-a",labl,
                                      "-a_srs",proj,
                                      "-te",str(geo[0]),str(geo[1]),str(geo[2]),str(geo[3]),
                                      name1,name2])

def assign_val(rec,df,label):

    rec = rec.copy()
    if rec['lon'].size == df['x'].size:
        rec[label] = np.nan
        for i in range(rec['lon'].size):
            dis = haversine.haversine_array(np.array(df['y'],dtype='float32'),
                                            np.array(df['x'],dtype='float32'),
                                            np.float32(rec.loc[i,'lat']),
                                            np.float32(rec.loc[i,'lon']))
            idx = np.argmin(dis)
            rec.loc[i,label] = df.loc[idx,label]

    return rec

def remove_outliers_qua(df,groby,label):
    
    df              = df.copy()
    df.loc[:,label] = df.groupby(groby)[label].transform(lambda x: x.where((x>x.quantile(0.01)) & (x<x.quantile(0.99)),np.nan))
    df.loc[:,label] = df.groupby(groby)[label].transform(pd.DataFrame.interpolate)
    df.loc[:,label] = df.groupby(groby)[label].fillna(method='ffill')
    df.loc[:,label] = df.groupby(groby)[label].fillna(method='bfill')
    # df.interpolate(inplace=True)

    return df

def remove_outliers_mad(df,groby,label):

    df              = df.copy()
    mask            = df.groupby(groby)[label].transform(is_outlier,thresh=3.5)
    df[label]       = df[label].where(~mask, np.nan)
    df.loc[:,label] = df.groupby(groby)[label].fillna(method='ffill')
    df.loc[:,label] = df.groupby(groby)[label].fillna(method='bfill')

    return df

def remove_outliers_yama(df,groby,label):

    df              = df.copy()
    df.loc[:,label] = df.groupby(groby)[label].transform(bank4flood)
    df.loc[:,label] = df.groupby(groby)[label].fillna(method='ffill')
    df.loc[:,label] = df.groupby(groby)[label].fillna(method='bfill')

    return df

def is_outlier(points, thresh=3.5):

    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """

    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh

if __name__ == '__main__':
    main(sys.argv[1:])
