#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: j.sosa@bristol.ac.uk / sosa.jeison@gmail.com

import os
import sys
import shutil
import getopt
import subprocess
import configparser
import pandas as pd
import numpy as np
import gdalutils
from osgeo import osr
from osgeo import gdal
from lfptools import shapefile
from lfptools.prepdata_utils import cy_d82d4
from lfptools.prepdata_utils import cy_rastermask
from lfptools.prepdata_utils import cy_directions_tau
from lfptools.prepdata_utils import cy_directions_esri
from lfptools.prepdata_utils import cy_rasterthreshold
from lfptools.prepdata_utils import calc_area


def prepdata(argv):
    """
    Info
    ----
    Delineate basins in the specified area

    Clip globally gridded distributed DEM, Hydrography and River Width data.

    Two parameters are required in this program: 1) Threshold and 2) Extent.
    Threshold to define river network in Km2. For example 5000 will delineate
    a river network with upslope area larger than 5000 Km2. Extent defines 
    clipped area. Four parameters should be indicated: xmin, ymin, xmax and ymax.


    Dependencies
    ------------
    TAUDEM : streamnet, gagewatershed
    GDAL : gdalwarp
    Pandas, Numpy, GDAL Python API, gdalutils, lfptools.prepdata_utils


    Inputs
    ------
    te : Clipping box specified as xmin,ymin,xmax,ymax (no space after commas!)
    out : Output folder
    dem : Any GDAL format (e.g. .tif, .vrt) containing DEM info
    acc : Any GDAL format (e.g. .tif, .vrt) containing accumulation info
    dir : Any GDAL format (e.g. .tif, .vrt) containing flow direction info
    thresh : Threshold to mask flow accumulation in KM**2
    streamnet : Calculate tree and coord files <yes/no>


    Outputs (If running at 30s)
    ---------------------------
    acc30.tif
    acc30_.tif
    area30.tif
    basins30.tif
    basins30d4.tif
    dem3.tif
    dir30.tif
    dir30tau.tif
    dir30tau_mask.tif
    dir30tau_maskd4.tif
    dir30taud4.tif
    net30.tif
    net30d4.tif
    out30.shp
    out30.tif
    out30d4.shp
    out30d4.tif
    stren_net30d4.out (folder)
    stren_net30d8.out (folder)
    stren_w30d4.tif
    stren_w30d8.tif
    strn_coord30d4.txt
    strn_coord30d8.txt
    strn_ord30d4.tif
    strn_ord30d8.tif
    strn_tree30d4.txt
    strn_tree30d8.txt
    """

    opts, args = getopt.getopt(argv, "i:")
    for o, a in opts:
        if o == "-i":
            inifile = a

    config = configparser.SafeConfigParser({'overwrite':False,'acc_area':False})
    config.read(inifile)

    te = np.float64(config.get('prepdata', 'te').split(','))
    out = str(config.get('prepdata', 'out'))
    _dem = str(config.get('prepdata', 'dem'))
    _acc = str(config.get('prepdata', 'acc'))
    _dir = str(config.get('prepdata', 'dir'))
    nproc = str(config.get('prepdata', 'nproc'))
    thresh = np.float64(config.get('prepdata', 'thresh'))
    streamnet = str(config.get('prepdata', 'streamnet'))
    overwrite = config.get('prepdata', 'overwrite').lower()=='True'.lower()
    acc_area = config.get('prepdata', 'acc_area').lower()=='True'.lower()

    # Defining extent
    xmin0 = te[0]
    ymin0 = te[1]
    xmax0 = te[2]
    ymax0 = te[3]

    # Getting resolution in degrees
    geo = gdalutils.get_geo(_acc)
    deg = round(geo[6], 4)

    # Checking resolution in arc-sec
    if deg == 0.0083:
        res = 30
    elif deg == 0.0008:
        res = 3

    # Creating an output folder
    try:
        os.makedirs(out)
    except OSError:
        if not os.path.isdir(out):
            raise

    # List of files generated
    dem3tif = out+'/dem3.tif'
    dir3tif = out+'/dir3.tif'

    dir3tau = out+'/dir3tau.tif'
    dir3taud4 = out+'/dir3taud4.tif'
    dir3tau_mask = out+'/dir3tau_mask.tif'
    dir3tau_maskd4 = out+'/dir3tau_maskd4.tif'
    dir30tif = out+'/dir30.tif'
    dir30tau = out+'/dir30tau.tif'
    dir30taud4 = out+'/dir30taud4.tif'
    dir30tau_mask = out+'/dir30tau_mask.tif'
    dir30tau_maskd4 = out+'/dir30tau_maskd4.tif'

    _acc3tif = out+'/acc3_.tif'
    acc3tif = out+'/acc3.tif'
    _acc30tif = out+'/acc30_.tif'
    acc30tif = out+'/acc30.tif'
    # _acc3tif is accumulation in grid cells, acc3tif is accumulation in area
    # If input acc is an area, override _acc3tif by the area version (then dont multiply by area)
    if acc_area:
        _acc3tif = acc3tif
        _acc30tif = acc30tif

    net3tif = out+'/net3.tif'
    net30tif = out+'/net30.tif'
    net3tifd4 = out+'/net3d4.tif'
    net30tifd4 = out+'/net30d4.tif'

    strn_ord3d8 = out+'/strn_ord3d8.tif'
    strn_tree3d8 = out+'/strn_tree3d8.txt'
    strn_coord3d8 = out+'/strn_coord3d8.txt'
    stren_net3d8 = out+'/stren_net3d8.out'
    stren_w3d8 = out+'/stren_w3d8.tif'
    strn_ord3d4 = out+'/strn_ord3d4.tif'
    strn_tree3d4 = out+'/strn_tree3d4.txt'
    strn_coord3d4 = out+'/strn_coord3d4.txt'
    stren_net3d4 = out+'/stren_net3d4.out'
    stren_w3d4 = out+'/stren_w3d4.tif'
    strn_ord30d8 = out+'/strn_ord30d8.tif'
    strn_tree30d8 = out+'/strn_tree30d8.txt'
    strn_coord30d8 = out+'/strn_coord30d8.txt'
    stren_net30d8 = out+'/stren_net30d8.out'
    stren_w30d8 = out+'/stren_w30d8.tif'
    strn_ord30d4 = out+'/strn_ord30d4.tif'
    strn_tree30d4 = out+'/strn_tree30d4.txt'
    strn_coord30d4 = out+'/strn_coord30d4.txt'
    stren_net30d4 = out+'/stren_net30d4.out'
    stren_w30d4 = out+'/stren_w30d4.tif'

    out3shp = out+'/out3.shp'
    out3shpd4 = out+'/out3d4.shp'
    out30shp = out+'/out30.shp'
    out30shpd4 = out+'/out30d4.shp'

    cat3tif = out+'/basins3.tif'
    cat3tifd4 = out+'/basins3d4.tif'
    cat30tif = out+'/basins30.tif'
    cat30tifd4 = out+'/basins30d4.tif'

    are3tif = out+'/area3.tif'
    are30tif = out+'/area30.tif'

    # Snap extent to match input tif grid cells
    geo = gdalutils.get_geo(_dem)
    # Geo has format [xmin, ymin, xmax, ymax, xn, yn, xres, yres, ....]
    xmin = geo[0] + np.floor((xmin0 - geo[0])/geo[6])*geo[6]
    ymin = geo[1] + np.floor((ymin0 - geo[1])/geo[7])*geo[7]
    xmax = geo[2] + np.floor((xmax0 - geo[2])/geo[6])*geo[6]
    ymax = geo[3] + np.floor((ymax0 - geo[3])/geo[7])*geo[7]

    # Clipping DEM .vrt files
    if not os.path.exists(dem3tif) or overwrite:
        print('clipping dem to region',xmin,ymin,xmax,ymax)
        subprocess.call(["gdalwarp", "-ot", "Float32", "-te", str(xmin), str(ymin), str(xmax),
                     str(ymax), "-overwrite", "-dstnodata", "-9999", "-co",'COMPRESS=DEFLATE',"-co", "BIGTIFF=YES", _dem, dem3tif])

    ########################################################################################
    # 3s resolution case
    #
    if res == 3:

        # Snap extent to match input tif grid cells
        geo = gdalutils.get_geo(_dir)
        # Geo has format [xmin, ymin, xmax, ymax, xn, yn, xres, yres, ....]
        xmin = geo[0] + np.floor((xmin0 - geo[0])/geo[6])*geo[6]
        ymin = geo[1] + np.floor((ymin0 - geo[1])/geo[7])*geo[7]
        xmax = geo[2] + np.floor((xmax0 - geo[2])/geo[6])*geo[6]
        ymax = geo[3] + np.floor((ymax0 - geo[3])/geo[7])*geo[7]
        print('clipping dir and acc fields to region',xmin,ymin,xmax,ymax)

        if not os.path.exists(dir3tif) or overwrite:
            subprocess.call(["gdalwarp", "-te", str(xmin), str(ymin), str(xmax),
                         str(ymax), "-overwrite", "-co", "BIGTIFF=YES","-co",'COMPRESS=DEFLATE', _dir, dir3tif])

        if not os.path.exists(_acc3tif) or overwrite:
            subprocess.call(["gdalwarp", "-te", str(xmin), str(ymin), str(xmax),
                         str(ymax), "-overwrite", "-co", "BIGTIFF=YES","-co",'COMPRESS=DEFLATE', _acc, _acc3tif])

        if not os.path.exists(dir3tau) or overwrite:
            print("converting directions into TAUDEM directions...")
            directions_tau(dir3tif, dir3tau)

        if not os.path.exists(are3tif) or overwrite:
            print("calculating area in extent...")
            calculate_area(dir3tau, are3tif)

        if not acc_area and (not os.path.exists(acc3tif) or overwrite):
            print("getting flow accumulation in km2...")
            multiply_rasters(_acc3tif, are3tif, acc3tif)

        if not os.path.exists(net3tif) or overwrite:
            print("thresholding accumulation to get river network...")
            rasterthreshold(acc3tif, thresh, 'Int16', net3tif)

        if not os.path.exists(dir3tau_mask) or overwrite:
            print("masking directions based on river network...")
            rastermask(dir3tau, net3tif, "Int16", dir3tau_mask)

        if not os.path.exists(out3shp) or overwrite:
            print("writing outlets and inland depressions in shapefile...")
            write_outlets(out3shp, dir3tau_mask)

        if not os.path.exists(cat3tif) or overwrite:
            print("writing basins file...")
            subprocess.call(["gagewatershed", "-p", dir3tau,
                         "-gw", cat3tif, "-o", out3shp])

        if streamnet == 'yes':
            # Streamnet fails if stren_net exists so remove first
            if os.path.exists(stren_net3d8) and overwrite: 
                shutil.rmtree(stren_net3d8)
            if not os.path.exists(stren_net3d8):
                # PFU: input -fel = dem3tif for correct slope in output streamnet
                subprocess.call(["mpiexec", "-n", nproc, "streamnet", "-fel", dem3tif, "-p", dir3tau, "-ad8", acc3tif, "-src", net3tif, "-ord",
                             strn_ord3d8, "-tree", strn_tree3d8, "-coord", strn_coord3d8, "-net", stren_net3d8, "-w", stren_w3d8, "-o", out3shp])

        if not os.path.exists(dir3tau_maskd4) or overwrite:
            print("creating D4 river network...")
            d82d4(dir3tau_mask, dir3tau_maskd4, net3tifd4)

        if not os.path.exists(out3shpd4) or overwrite:
            print("writing D4 outlets and inland depression in shapefile")
            write_outlets(out3shpd4, dir3tau_maskd4)

        if not os.path.exists(dir3taud4) or overwrite:
            print("create flow directions map D4...")
            create_dir_d4(dir3taud4, dir3tau, dir3tau_maskd4)

        if not os.path.exists(cat3tifd4) or overwrite:
            print("writing basins file D4...")
            subprocess.call(["gagewatershed", "-p", dir3taud4,
                         "-gw", cat3tifd4, "-o", out3shpd4])

        if streamnet == 'yes':
            # Streamnet fails if stren_net exists so remove first
            if os.path.exists(stren_net3d4) and overwrite: 
                shutil.rmtree(stren_net3d4)
            if not os.path.exists(stren_net3d4):
                # PFU: input -fel = dem3tif for correct slope in output streamnet
                subprocess.call(["mpiexec", "-n", nproc, "streamnet", "-fel", dem3tif, "-p", dir3tau_maskd4, "-ad8", acc3tif, "-src", net3tifd4,
                             "-ord", strn_ord3d4, "-tree", strn_tree3d4, "-coord", strn_coord3d4, "-net", stren_net3d4, "-w", stren_w3d4, "-o", out3shpd4])


    ########################################################################################
    # 30s resolution case
    #
    elif res == 30:
        # Snap extent to match input tif grid cells
        geo = gdalutils.get_geo(_dir)
        # Geo has format [xmin, ymin, xmax, ymax, xn, yn, xres, yres, ....]
        xmin = geo[0] + np.floor((xmin0 - geo[0])/geo[6])*geo[6]
        ymin = geo[1] + np.floor((ymin0 - geo[1])/geo[7])*geo[7]
        xmax = geo[2] + np.floor((xmax0 - geo[2])/geo[6])*geo[6]
        ymax = geo[3] + np.floor((ymax0 - geo[3])/geo[7])*geo[7]
        if not os.path.exists(dir30tif) or overwrite:
            subprocess.call(["gdalwarp", "-te", str(xmin), str(ymin),
                         str(xmax), str(ymax), "-overwrite", _dir, dir30tif])

        if not os.path.exists(_acc30tif) or overwrite:
            subprocess.call(["gdalwarp", "-te", str(xmin), str(ymin), str(xmax),
                         str(ymax), "-overwrite", "-co", "BIGTIFF=YES", _acc, _acc30tif])

        if not os.path.exists(dir30tau) or overwrite:
            print("converting directions into TAUDEM directions...")
            directions_tau(dir30tif, dir30tau)

        if not os.path.exists(are30tif) or overwrite:
            print("calculating area in extent...")
            calculate_area(dir30tau, are30tif)

        if not acc_area and (not os.path.exists(acc30tiff) or overwrite):
            print("getting flow accumulation in km2...")
            multiply_rasters(_acc30tif, are30tif, acc30tif)

        if not os.path.exists(net30tif) or overwrite:
            print("thresholding accumulation to get river network...")
            rasterthreshold(acc30tif, thresh, 'Int16', net30tif)

        if not os.path.exists(dir30tau_mask) or overwrite:
            print("masking directions based on river network...")
            rastermask(dir30tau, net30tif, "Int16", dir30tau_mask)

        if not os.path.exists(out30shp) or overwrite:
            print("writing outlets and inland depressions in shapefile...")
            write_outlets(out30shp, dir30tau_mask)

        if not os.path.exists(cat30tif) or overwrite:
            print("writing basins file...")
            subprocess.call(["gagewatershed", "-p", dir30tau,
                         "-gw", cat30tif, "-o", out30shp])

        if streamnet == 'yes':
            # Streamnet fails if stren_net exists so remove first
            if os.path.exists(stren_net30d8) and overwrite: 
                shutil.rmtree(stren_net30d8)
            if not os.path.exists(stren_net30d8):
                # PFU: input -fel should be dem for correct slope in output stremnet
                # BUT we dont have a dem file at 30s
                subprocess.call(["mpiexec", "-n", nproc, "streamnet", "-fel", net30tif, "-p", dir30tau, "-ad8", acc30tif, "-src", net30tif, "-ord",
                             strn_ord30d8, "-tree", strn_tree30d8, "-coord", strn_coord30d8, "-net", stren_net30d8, "-w", stren_w30d8, "-o", out30shp])

        if not os.path.exists(dir30tau_maskd4) or overwrite:
            print("creating D4 river network...")
            d82d4(dir30tau_mask, dir30tau_maskd4, net30tifd4)

        if not os.path.exists(out30shpd4) or overwrite:
            print("writing D4 outlets and inland depression in shapefile...")
            write_outlets(out30shpd4, dir30tau_maskd4)

        if not os.path.exists(dir30taud4) or overwrite:
            print("create flow directions map D4...")
            create_dir_d4(dir30taud4, dir30tau, dir30tau_maskd4)

        if not os.path.exists(cat30tifd4) or overwrite:
            print("writing basins file D4...")
            subprocess.call(["gagewatershed", "-p", dir30taud4,
                         "-gw", cat30tifd4, "-o", out30shpd4])

        if streamnet == 'yes':
            # Streamnet fails if stren_net exists so remove first
            if os.path.exists(stren_net30d4) and overwrite: 
                shutil.rmtree(stren_net30d4)
            if not os.path.exists(stren_net30d4):
                # PFU: input -fel should be dem for correct slope in output stremnet
                # BUT we dont have a dem file at 30s
                subprocess.call(["mpiexec", "-n", nproc, "streamnet", "-fel", net30tifd4, "-p", dir30tau_maskd4, "-ad8", acc30tif, "-src", net30tifd4,
                             "-ord", strn_ord30d4, "-tree", strn_tree30d4, "-coord", strn_coord30d4, "-net", stren_net30d4, "-w", stren_w30d4, "-o", out30shpd4])


def directions_tau(inputrast, outputrast):
    """
    Function to use in Shell to change convetion from a DIR file
    HydroSHEDS uses ESRI convention 128,64,32,.. this script
    changes these numbers to TauDEM convention, 1,2,3...
    """

    nodata = -32768
    data = gdalutils.get_data(inputrast)
    datageo = gdalutils.get_geo(inputrast)
    datatau = cy_directions_tau(np.int16(data), np.int16(nodata))
    gdalutils.write_raster(np.float64(
        datatau), outputrast, datageo, "Int16", nodata)

def directions_esri(inputrast, outputrast):
    """
    Function to change convetion from a DIR file
    This script changes these numbers from TauDEM convention, 1,2,3...
	to ESRI convention 128,64,32,.. 
    """

    nodata = -32768
    data = gdalutils.get_data(inputrast)
    datageo = gdalutils.get_geo(inputrast)
    data_esri = cy_directions_esri(np.int16(data), np.int16(nodata))
    gdalutils.write_raster(np.int16(data_esri), outputrast, datageo, "Int16", nodata)


def rasterthreshold(file, thres, fmt, outp):
    """
    Output a raster based on a threshold (larger-equal-than)
    """

    nodata = -1
    filedata = gdalutils.get_data(file)
    filegeo = gdalutils.get_geo(file)
    data = cy_rasterthreshold(np.float64(
        filedata), np.float64(thres), np.float64(nodata))
    gdalutils.write_raster(np.float64(data), outp, filegeo, fmt, nodata)


def rastermask(file, mask, fmt, outp):
    """
    Mask input array following a defined mask (1,0)
    """

    nodata = -32768
    filedata = gdalutils.get_data(file)
    maskdata = gdalutils.get_data(mask)
    filegeo = gdalutils.get_geo(file)
    data = cy_rastermask(np.float64(filedata), np.int16(maskdata))
    gdalutils.write_raster(np.float64(data), outp, filegeo, fmt, nodata)


def mosaic_region(inputpath, xmin, ymin, xmax, ymax, outputfile):

    mytiles = listdir(inputpath, '.tif')
    myoutput = open(outputfile, 'w')

    for i in range(len(mytiles)):

        fname = os.path.basename(mytiles[i])

        if fname[0] == 'n':
            lat1 = np.float64(fname[1:3])
        else:
            lat1 = -np.float64(fname[1:3])
        lat2 = lat1+5

        if fname[3] == 'e':
            lon1 = np.float64(fname[4:7])
        else:
            lon1 = -np.float64(fname[4:7])
        lon2 = lon1+5

        if (lat2 <= ymax) & (lat1 >= ymin) & (lon2 <= xmax) & (lon1 >= xmin):
            myoutput.write(mytiles[i]+'\n')

    myoutput.close()


def d82d4(filedir, filediro, fileneto):
    """
    Returns direction and river netwrok maps in D4
    """

    nodata = -32768.
    dirdata = gdalutils.get_data(filedir)
    dirgeo = gdalutils.get_geo(filedir)
    data, net = cy_d82d4(np.int16(dirdata), np.int16(nodata))
    gdalutils.write_raster(np.int16(data), filediro, dirgeo, "Int16", nodata)
    gdalutils.write_raster(np.int16(net), fileneto, dirgeo, "Int16", nodata)


def listdir(path, ext):
    """ Returns a list of files in a folder based on a specified extension """

    mylist = []
    for root, dirs, files in os.walk(path):
        for filename in files:
            if filename.endswith(ext):
                mylist.append(os.path.join(root, filename))
    return mylist


def write_list_files(inputpath, ext, outputfile):

    mytiles = listdir(inputpath, ext)
    myoutput = open(outputfile, 'w')

    for i in range(len(mytiles)):
        myoutput.write(mytiles[i]+'\n')
    myoutput.close()


def write_outlets(outshp, dirtif_mask):

    proj = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

    dat = gdalutils.get_data(dirtif_mask)
    geo = gdalutils.get_geo(dirtif_mask)
    rows, cols = np.where(dat > 0)

    x = []
    y = []
    for row, col in zip(rows, cols):
        A = find_neighbours(dat, row, col)
        if np.any(A < 0):
            x.append(geo[8][col])
            y.append(geo[9][row])

    # Initiate shapefile
    w = shapefile.Writer(shapefile.POINT)
    w.field('x')
    w.field('y')
    w.field('id')

    # Write coordinate points in shapefile
    for i in range(len(x)):
        w.point(x[i], y[i])
        w.record(x[i], y[i], i)
    w.save(outshp)
    fname = os.path.dirname(outshp)+'/' + \
        os.path.basename(outshp).split('.')[0] + '.prj'
    prj = open(fname, "w")
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj)
    prj.write(srs.ExportToWkt())
    prj.close()

    typ = "Byte"
    fmt = "GTiff"
    nodata = 0
    name1 = os.path.dirname(outshp)+'/' + \
        os.path.basename(outshp).split('.')[0] + '.shp'
    name2 = os.path.dirname(outshp)+'/' + \
        os.path.basename(outshp).split('.')[0] + '.tif'
    subprocess.call(["gdal_rasterize", "-a_nodata", str(nodata), "-ot", typ, "-of", fmt, "-tr", str(geo[6]), str(geo[7]),
                     "-burn", "1", "-a_srs", proj, "-te", str(geo[0]), str(geo[1]), str(geo[2]), str(geo[3]), name1, name2])


def find_neighbours(dat, row, col):

    nei = []
    try:
        nei.append(dat[row, col+1])
    except IndexError:
        nei.append(9999)
    try:
        nei.append(dat[row-1, col+1])
    except IndexError:
        nei.append(9999)
    try:
        nei.append(dat[row-1, col])
    except IndexError:
        nei.append(9999)
    try:
        nei.append(dat[row-1, col-1])
    except IndexError:
        nei.append(9999)
    try:
        nei.append(dat[row, col-1])
    except IndexError:
        nei.append(9999)
    try:
        nei.append(dat[row+1, col-1])
    except IndexError:
        nei.append(9999)
    try:
        nei.append(dat[row+1, col])
    except IndexError:
        nei.append(9999)
    try:
        nei.append(dat[row+1, col+1])
    except IndexError:
        nei.append(9999)

    return np.array(nei)


def create_dir_d4(dirtaud4, dirtaud8, dirtau_maskd4):

    dat1 = gdalutils.get_data(dirtaud8)
    dat2 = gdalutils.get_data(dirtau_maskd4)
    geo = gdalutils.get_geo(dirtaud8)

    A = np.where(dat2 > 0)
    dat1[A] = dat2[A]

    gdalutils.write_raster(dat1, dirtaud4, geo, "Int16", -32768)


def read_tree_taudem(treef):
    df = pd.read_csv(treef, sep='\t', names=[
                     '0', 'link_no', 'start_pnt', 'end_pnt', 'frst_ds', 'frst_us', 'scnd_us', 'strahler', 'mon_pnt', 'shreve'])
    df.drop(['0'], axis=1, inplace=True)
    return df


def read_coord_taudem(coordf):
    df = pd.read_csv(coordf, sep='\t', names=[
                     '0', 'lon', 'lat', 'distance', 'elev', 'contr_area'])
    df.drop(['0'], axis=1, inplace=True)
    return df


def calculate_area(filename, output):

    geo = gdalutils.get_geo(filename)
    nx = np.int32(geo[4])
    ny = np.int32(geo[5])
    resx = np.float32(geo[6])
    resy = np.float32(geo[7])
    x = np.float32(geo[8])
    y = np.float32(geo[9])
    dat = calc_area(nx, ny, resx, resy, x, y)
    gdalutils.write_raster(np.array(dat), output, geo, "Float32", -9999)


def multiply_rasters(rast1, rast2, out):

    dat1 = gdalutils.get_data(rast1)
    dat2 = gdalutils.get_data(rast2)
    geo1 = gdalutils.get_geo(rast1)
    geo2 = gdalutils.get_geo(rast2)

    dat_masked1 = np.ma.masked_where(dat1 == geo1[11], dat1)
    dat_masked2 = np.ma.masked_where(dat2 == geo2[11], dat2)

    res = dat_masked1 * dat_masked2
    res.set_fill_value(-9999)

    gdalutils.write_raster(res.filled(), out, geo1, "Float32", -9999)


if __name__ == '__main__':
    prepdata(sys.argv[1:])
