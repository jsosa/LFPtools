#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# date: 22/nov/2017
# mail: j.sosa@bristol.ac.uk / sosa.jeison@gmail.com

import os
import sys
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
from lfptools.prepdata_utils import cy_rasterthreshold
from lfptools.prepdata_utils import calc_area


def prepdata(argv):

    """
    Program to extract Digital Elevation Model, Flow direction map and River width
    from Yamazaki's datasets. This program clip flow directions based on a specified
    box, then flow accumulation is calcualted by using TauDEM. A threshold is applied
    to define a river network mask

    Dependencies
    ------------

    GDAL : gdalwarp
    TauDEM : streamnet

    Parameters
    ----------

    <predata.cfg> file containing this text:

    [prepdata]

    te        = clipping box specified as xmin,ymin,xmax,ymax (no space after comma!)
    out       = output folder
    dem       = Any GDAL format (e.g. .tif, .vrt) containing DEM info
    acc       = Any GDAL format (e.g. .tif, .vrt) containing accumulation info
    dir       = Any GDAL format (e.g. .tif, .vrt) containing flow direction info
    wth       = Any GDAL format (e.g. .tif, .vrt) containing river width info
    thresh    = threshold to mask flow accumulation in KM**2
    streamnet = calculate tree and coord files <yes/no>

    Usage
    -----
    $ lfp-prepdata -i params.cfg
    """

    opts, args = getopt.getopt(argv,"i:")
    for o, a in opts:
        if o == "-i": inifile  = a

    config = configparser.SafeConfigParser()
    config.read(inifile)

    te        = np.float64(config.get('prepdata','te').split(','))
    out       = str(config.get('prepdata','out'))
    _dem      = str(config.get('prepdata','dem'))
    _acc      = str(config.get('prepdata','acc'))
    _dir      = str(config.get('prepdata','dir'))
    _wth      = str(config.get('prepdata','wth'))
    nproc     = str(config.get('prepdata','nproc'))
    thresh    = np.float64(config.get('prepdata','thresh'))
    streamnet = str(config.get('prepdata','streamnet'))
    
    # Defining extent
    xmin = te[0]
    ymin = te[1]
    xmax = te[2]
    ymax = te[3]

    # Getting resolution in degrees
    geo = gdalutils.get_geo(_acc)
    deg = round(geo[6],4)

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
    dem3tif         = './'+out+'/dem3.tif'
    wth3tif         = './'+out+'/wth3.tif'
    dir3tif         = './'+out+'/dir3.tif'

    dir3tau         = './'+out+'/dir3tau.tif'
    dir3tau_mask    = './'+out+'/dir3tau_mask.tif'
    dir3tau_maskd4  = './'+out+'/dir3tau_maskd4.tif'
    dir30tif        = './'+out+'/dir30.tif'
    dir30tau        = './'+out+'/dir30tau.tif'
    dir30taud4      = './'+out+'/dir30taud4.tif'
    dir30tau_mask   = './'+out+'/dir30tau_mask.tif'
    dir30tau_maskd4 = './'+out+'/dir30tau_maskd4.tif'
    
    _acc3tif        = './'+out+'/acc3_.tif'
    acc3tif         = './'+out+'/acc3.tif'
    _acc30tif       = './'+out+'/acc30_.tif'
    acc30tif        = './'+out+'/acc30.tif'


    net3tif         = './'+out+'/net3.tif'
    net30tif        = './'+out+'/net30.tif'
    net3tifd4       = './'+out+'/net3d4.tif'
    net30tifd4      = './'+out+'/net30d4.tif'

    strn_ord3d8     = './'+out+'/strn_ord3d8.tif'
    strn_tree3d8    = './'+out+'/strn_tree3d8.txt'
    strn_coord3d8   = './'+out+'/strn_coord3d8.txt'
    stren_net3d8    = './'+out+'/stren_net3d8.out'
    stren_w3d8      = './'+out+'/stren_w3d8.tif'
    strn_ord3d4     = './'+out+'/strn_ord3d4.tif'
    strn_tree3d4    = './'+out+'/strn_tree3d4.txt'
    strn_coord3d4   = './'+out+'/strn_coord3d4.txt'
    stren_net3d4    = './'+out+'/stren_net3d4.out'
    stren_w3d4      = './'+out+'/stren_w3d4.tif'
    strn_ord30d8    = './'+out+'/strn_ord30d8.tif'
    strn_tree30d8   = './'+out+'/strn_tree30d8.txt'
    strn_coord30d8  = './'+out+'/strn_coord30d8.txt'
    stren_net30d8   = './'+out+'/stren_net30d8.out'
    stren_w30d8     = './'+out+'/stren_w30d8.tif'
    strn_ord30d4    = './'+out+'/strn_ord30d4.tif'
    strn_tree30d4   = './'+out+'/strn_tree30d4.txt'
    strn_coord30d4  = './'+out+'/strn_coord30d4.txt'
    stren_net30d4   = './'+out+'/stren_net30d4.out'
    stren_w30d4     = './'+out+'/stren_w30d4.tif'

    out3shp         = './'+out+'/out3.shp'
    out3shpd4       = './'+out+'/out3d4.shp'
    out30shp        = './'+out+'/out30.shp'
    out30shpd4      = './'+out+'/out30d4.shp'

    cat3tif         = './'+out+'/basins3.tif'
    cat30tif        = './'+out+'/basins30.tif'
    cat30tifd4      = './'+out+'/basins30d4.tif'

    are3tif         = './'+out+'/area3.tif'
    are30tif        = './'+out+'/area30.tif'

    # Clipping DEM and river width .vrt files
    subprocess.call(["gdalwarp","-ot","Float32","-te",str(xmin),str(ymin),str(xmax),str(ymax),"-overwrite","-dstnodata","-9999","-co","BIGTIFF=YES",_dem,dem3tif])
    subprocess.call(["gdalwarp","-ot","Float32","-te",str(xmin),str(ymin),str(xmax),str(ymax),"-overwrite","-dstnodata","-9999","-co","BIGTIFF=YES",_wth,wth3tif])

    if res == 3:

        subprocess.call(["gdalwarp","-te",str(xmin),str(ymin),str(xmax),str(ymax),"-overwrite","-co","BIGTIFF=YES",_dir,dir3tif])

        subprocess.call(["gdalwarp","-te",str(xmin),str(ymin),str(xmax),str(ymax),"-overwrite","-co","BIGTIFF=YES",_acc,_acc3tif])

        print("converting directions into TAUDEM directions...")
        directions_tau(dir3tif,dir3tau)

        print("calculating area in extent...")
        calculate_area(dir3tau,are3tif)

        print("getting flow accumulation in km2...")
        multiply_rasters(_acc3tif,are3tif,acc3tif)

        print("thresholding accumulation to get river network...")
        rasterthreshold(acc3tif,thresh,'Int16',net3tif)

        print("masking directions based on river network...")
        rastermask(dir3tau,net3tif,"Int16",dir3tau_mask)

        print("writing outlets and inland depressions in shapefile...")
        write_outlets(out3shp,dir3tau_mask)

        print("writing basins file...")
        subprocess.call(["gagewatershed","-p",dir3tau,"-gw",cat3tif,"-o",out3shp])

        if streamnet == 'yes':
            subprocess.call(["mpiexec","-n",nproc,"streamnet","-fel",net3tif,"-p",dir3tau,"-ad8",acc3tif,"-src",net3tif,"-ord",strn_ord3d8,"-tree",strn_tree3d8,"-coord",strn_coord3d8,"-net",stren_net3d8,"-w",stren_w3d8,"-o",out3shp])

        print("creating D4 river network...")
        d82d4(dir3tau_mask,dir3tau_maskd4,net3tifd4)

        print("writing D4 outlets and inland depression in shapefile")
        write_outlets(out3shpd4,dir3tau_maskd4)

        print("create flow directions map D4...")
        create_dir_d4(dir3taud4,dir3tau,dir3tau_maskd4)

        print("writing basins file D4...")
        subprocess.call(["gagewatershed","-p",dir3taud4,"-gw",cat3tifd4,"-o",out3shpd4])

        if streamnet == 'yes':
            subprocess.call(["mpiexec","-n",nproc,"streamnet","-fel",net3tifd4,"-p",dir3tau_maskd4,"-ad8",acc3tif,"-src",net3tifd4,"-ord",strn_ord3d4,"-tree",strn_tree3d4,"-coord",strn_coord3d4,"-net",stren_net3d4,"-w",stren_w3d4,"-o",out3shpd4])

    elif res == 30:

        subprocess.call(["gdalwarp","-te",str(xmin),str(ymin),str(xmax),str(ymax),"-overwrite",_dir,dir30tif])

        subprocess.call(["gdalwarp","-te",str(xmin),str(ymin),str(xmax),str(ymax),"-overwrite","-co","BIGTIFF=YES",_acc,_acc30tif])

        print("converting directions into TAUDEM directions...")
        directions_tau(dir30tif,dir30tau)

        print("calculating area in extent...")
        calculate_area(dir30tau,are30tif)

        print("getting flow accumulation in km2...")
        multiply_rasters(_acc30tif,are30tif,acc30tif)

        print("thresholding accumulation to get river network...")
        rasterthreshold(acc30tif,thresh,'Int16',net30tif)

        print("masking directions based on river network...")
        rastermask(dir30tau,net30tif,"Int16",dir30tau_mask)

        print("writing outlets and inland depressions in shapefile...")
        write_outlets(out30shp,dir30tau_mask)

        print("writing basins file...")
        subprocess.call(["gagewatershed","-p",dir30tau,"-gw",cat30tif,"-o",out30shp])

        if streamnet == 'yes':
            subprocess.call(["mpiexec","-n",nproc,"streamnet","-fel",net30tif,"-p",dir30tau,"-ad8",acc30tif,"-src",net30tif,"-ord",strn_ord30d8,"-tree",strn_tree30d8,"-coord",strn_coord30d8,"-net",stren_net30d8,"-w",stren_w30d8,"-o",out30shp])

        print("creating D4 river network...")
        d82d4(dir30tau_mask,dir30tau_maskd4,net30tifd4)

        print("writing D4 outlets and inland depression in shapefile...")
        write_outlets(out30shpd4,dir30tau_maskd4)

        print("create flow directions map D4...")
        create_dir_d4(dir30taud4,dir30tau,dir30tau_maskd4)

        print("writing basins file D4...")
        subprocess.call(["gagewatershed","-p",dir30taud4,"-gw",cat30tifd4,"-o",out30shpd4])

        if streamnet == 'yes':
            subprocess.call(["mpiexec","-n",nproc,"streamnet","-fel",net30tifd4,"-p",dir30tau_maskd4,"-ad8",acc30tif,"-src",net30tifd4,"-ord",strn_ord30d4,"-tree",strn_tree30d4,"-coord",strn_coord30d4,"-net",stren_net30d4,"-w",stren_w30d4,"-o",out30shpd4])


def directions_tau(inputrast,outputrast):

    """
    Function to use in Shell to change convetion from a DIR file
    HydroSHEDS uses ESRI convention 128,64,32,.. this script
    changes these numbers to TauDEM convention, 1,2,3...
    """

    nodata  = -32768
    data    = gdalutils.get_data(inputrast)
    datageo = gdalutils.get_geo(inputrast)
    datatau = cy_directions_tau(np.int16(data),np.int16(nodata))
    gdalutils.write_raster(np.float64(datatau),outputrast,datageo,"Int16",nodata)


def rasterthreshold(file,thres,fmt,outp):

    """
    Output a raster based on a threshold (larger-equal-than)
    """

    nodata   = -1
    filedata = gdalutils.get_data(file)
    filegeo  = gdalutils.get_geo(file)
    data     = cy_rasterthreshold(np.float64(filedata), np.float64(thres), np.float64(nodata))
    gdalutils.write_raster(np.float64(data),outp,filegeo,fmt,nodata)


def rastermask(file,mask,fmt,outp):

    """
    Mask input array following a defined mask (1,0)
    """

    nodata = -32768
    filedata = gdalutils.get_data(file)
    maskdata = gdalutils.get_data(mask)
    filegeo  = gdalutils.get_geo(file)
    data     = cy_rastermask(np.float64(filedata),np.int16(maskdata))
    gdalutils.write_raster(np.float64(data),outp,filegeo,fmt,nodata)


def mosaic_region(inputpath,xmin,ymin,xmax,ymax,outputfile):

    mytiles  = listdir(inputpath,'.tif')
    myoutput = open(outputfile,'w')

    for i in range(len(mytiles)):

        fname = os.path.basename(mytiles[i])

        if fname[0] == 'n': lat1 = np.float64(fname[1:3])
        else: lat1 = -np.float64(fname[1:3])
        lat2  = lat1+5

        if fname[3] == 'e': lon1 = np.float64(fname[4:7])
        else: lon1 = -np.float64(fname[4:7])
        lon2  = lon1+5

        if (lat2<=ymax) & (lat1>=ymin) & (lon2<=xmax) & (lon1>=xmin):
            myoutput.write(mytiles[i]+'\n')

    myoutput.close()


def d82d4(filedir,filediro,fileneto):

    """
    Returns direction and river netwrok maps in D4
    """

    nodata   = -32768.
    dirdata  = gdalutils.get_data(filedir)
    dirgeo   = gdalutils.get_geo(filedir)
    data,net = cy_d82d4(np.int16(dirdata), np.int16(nodata))
    gdalutils.write_raster(np.int16(data),filediro,dirgeo,"Int16",nodata)
    gdalutils.write_raster(np.int16(net),fileneto,dirgeo,"Int16",nodata)


def listdir(path,ext):

    """ Returns a list of files in a folder based on a specified extension """

    mylist = []
    for root, dirs, files in os.walk(path):
        for filename in files:
            if filename.endswith(ext):
                mylist.append(os.path.join(root, filename))
    return mylist


def write_list_files(inputpath,ext,outputfile):

    mytiles  = listdir(inputpath,ext)
    myoutput = open(outputfile,'w')

    for i in range(len(mytiles)):
        myoutput.write(mytiles[i]+'\n')
    myoutput.close()


def write_outlets(outshp,dirtif_mask):

    proj = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

    dat = gdalutils.get_data(dirtif_mask)
    geo = gdalutils.get_geo(dirtif_mask)
    rows,cols = np.where(dat>0)

    x = []
    y = []
    for row,col in zip(rows,cols):
        A=find_neighbours(dat,row,col)
        if np.any(A<0):
            x.append(geo[8][col])
            y.append(geo[9][row])

    # Initiate shapefile
    w = shapefile.Writer(shapefile.POINT)
    w.field('x')
    w.field('y')
    w.field('id')

    # Write coordinate points in shapefile
    for i in range(len(x)):
        w.point(x[i],y[i])
        w.record(x[i],y[i],i)
    w.save(outshp)
    fname = os.path.dirname(outshp)+'/'+os.path.basename(outshp).split('.')[0] + '.prj'
    prj = open(fname, "w")
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj)
    prj.write(srs.ExportToWkt())
    prj.close()

    typ    = "Byte"
    fmt    = "GTiff"
    nodata = 0
    name1  = os.path.dirname(outshp)+'/'+os.path.basename(outshp).split('.')[0] + '.shp'
    name2  = os.path.dirname(outshp)+'/'+os.path.basename(outshp).split('.')[0] + '.tif'
    subprocess.call(["gdal_rasterize","-a_nodata",str(nodata),"-ot",typ,"-of",fmt,"-tr",str(geo[6]),str(geo[7]),"-burn","1","-a_srs",proj,"-te",str(geo[0]),str(geo[1]),str(geo[2]),str(geo[3]),name1,name2])


def find_neighbours(dat,row,col):

    nei = []
    try:
        nei.append(dat[row,col+1])
    except IndexError:
        nei.append(9999)
    try:
        nei.append(dat[row-1,col+1])
    except IndexError:
        nei.append(9999)
    try:
        nei.append(dat[row-1,col])
    except IndexError:
        nei.append(9999)
    try:
        nei.append(dat[row-1,col-1])
    except IndexError:
        nei.append(9999)
    try:
        nei.append(dat[row,col-1])
    except IndexError:
        nei.append(9999)
    try:
        nei.append(dat[row+1,col-1])
    except IndexError:
        nei.append(9999)
    try:
        nei.append(dat[row+1,col])
    except IndexError:
        nei.append(9999)
    try:
        nei.append(dat[row+1,col+1])
    except IndexError:
        nei.append(9999)

    return np.array(nei)


def create_dir_d4(dirtaud4,dirtaud8,dirtau_maskd4):

    dat1 = gdalutils.get_data(dirtaud8)
    dat2 = gdalutils.get_data(dirtau_maskd4)
    geo  = gdalutils.get_geo(dirtaud8)

    A = np.where(dat2>0)
    dat1[A] = dat2[A]

    gdalutils.write_raster(dat1,dirtaud4,geo,"Int16",-32768)


def read_tree_taudem(treef):
    df = pd.read_csv(treef, sep='\t', names=['0','link_no','start_pnt','end_pnt','frst_ds','frst_us','scnd_us','strahler','mon_pnt','shreve'])
    df.drop(['0'], axis=1, inplace=True)
    return df


def read_coord_taudem(coordf):
    df = pd.read_csv(coordf, sep='\t', names=['0','lon','lat','distance','elev','contr_area'])
    df.drop(['0'], axis=1, inplace=True)
    return df


def calculate_area(filename,output):

    geo  = gdalutils.get_geo(filename)
    nx   = np.int16(geo[4])
    ny   = np.int16(geo[5])
    resx = np.float32(geo[6])
    resy = np.float32(geo[7])
    x    = np.float32(geo[8])
    y    = np.float32(geo[9])
    dat  = calc_area(nx,ny,resx,resy,x,y)
    gdalutils.write_raster(np.array(dat),output,geo,"Float32",-9999)


def multiply_rasters(rast1,rast2,out):

    dat1 = gdalutils.get_data(rast1)
    dat2 = gdalutils.get_data(rast2)
    geo1 = gdalutils.get_geo(rast1)
    geo2 = gdalutils.get_geo(rast2)
    
    dat_masked1 = np.ma.masked_where(dat1==geo1[11],dat1)
    dat_masked2 = np.ma.masked_where(dat2==geo2[11],dat2)

    res = dat_masked1 * dat_masked2
    res.set_fill_value(-9999)

    gdalutils.write_raster(res.filled(),out,geo1,"Float32",-9999)

if __name__ == '__main__':
    prepdata(sys.argv[1:])
