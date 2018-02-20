#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# date: 22/nov/2017
# mail: j.sosa@bristol.ac.uk / sosa.jeison@gmail.com

import os
import sys
import getopt
import subprocess
import shapefile
import pandas as pd
import numpy as np
import gdalutils
from osgeo import osr
from osgeo import gdal
from prepdata_utils import cy_d82d4
from prepdata_utils import cy_rastermask
from prepdata_utils import cy_directions_tau
from prepdata_utils import cy_rasterthreshold
from prepdata_utils import calc_area


def prepdata(argv):

    """
    Program to extract Digital Elevation Model, Flow direction map and River width
    from Yamazaki's datasets. This program clip flow directions based on a specified
    box, then flow accumulation is calcualted by using TauDEM. A threshold is applied
    to define a river network mask

    Dependencies
    ------------

    GDAL : gdalbuildvrt, gdalwarp
    TauDEM : aread8, streamnet

    Parameters
    ----------

    --te : clipping box specified as xmin,ymin,xmax,ymax (no space after comma!)
    --out : output folder
    --srcdir : flow directions source <yamazaki, hydrosheds>
    --resdir : flow directions resolution <3,30>
    --thresh : threshold to mask flow accumulation
    --srcdem : DEM source <yamazaki,hydrosheds>
    --streamnet : calculate tree and coord files <yes/no>

    Usage
    -----
    $ python prepdata.py --te -30,30,45,72 --out eu_yamazaki --srcdir yamazaki --srcdem yamazaki --resdir 3 --thresh 100000 --streamnet no

    """

    opts, args = getopt.getopt(argv,"a:b:c:d:",["te=","out=","resdir=","srcdir=","srcdem=","thresh=","streamnet="])
    for o, a in opts:
        if o == "--te"        : te        = np.float64(a.split(','))
        if o == "--out"       : outf      = str(a)
        if o == "--srcdir"    : srcdir    = str(a)
        if o == "--resdir"    : res       = np.int(a)
        if o == "--thresh"    : thresh    = np.float64(a)
        if o == "--srcdem"    : srcdem    = str(a)
        if o == "--streamnet" : streamnet = str(a)

    xmin = te[0]
    ymin = te[1]
    xmax = te[2]
    ymax = te[3]

    try:
        os.makedirs(outf)
    except OSError:
        if not os.path.isdir(outf):
            raise

    # source files
    if srcdem == "yamazaki":
        dem3path = '/scratch/hydro3/js16170/data_processed/yamazaki/dem3/'
    elif srcdem == "hydrosheds":
        dem3path = '/scratch/hydro3/js16170/data_processed/hydrosheds/dem3/'

    if (srcdir == "yamazaki") & (res == 3):
        dir3path = '/scratch/hydro3/js16170/data_processed/yamazaki/dir3/'
    elif (srcdir == "hydrosheds") & (res == 3):
        dir3path = '/scratch/hydro3/js16170/data_processed/hydrosheds/dir3/'

    if (srcdir == "yamazaki") & (res == 30):
        sys.exit("ERROR: There isn't Yamazaki's 30s Flow direction data")
    elif (srcdir == "hydrosheds") & (res == 30):
        dir30path = '/scratch/hydro3/js16170/data_processed/hydrosheds/dir30/'
        dir30merged = dir30path + "merged.tif"

    wth3path = '/scratch/hydro3/js16170/data_processed/yamazaki/wth3/'

    # list of files generated
    dem3list        = './'+outf+'/dem3.txt'
    dem3vrt         = './'+outf+'/dem3.vrt'
    dem3tif         = './'+outf+'/dem3.tif'
    wth3list        = './'+outf+'/wth3.txt'
    wth3vrt         = './'+outf+'/wth3.vrt'
    wth3tif         = './'+outf+'/wth3.tif'
    dir3list        = './'+outf+'/dir3.txt'
    dir3vrt         = './'+outf+'/dir3.vrt'
    dir3tif         = './'+outf+'/dir3.tif'
    dir3tau         = './'+outf+'/dir3tau.tif'
    dir3tau_mask    = './'+outf+'/dir3tau_mask.tif'
    dir3tau_maskd4  = './'+outf+'/dir3tau_maskd4.tif'
    dir30list       = './'+outf+'/dir30.txt'
    dir30vrt        = './'+outf+'/dir30.vrt'
    dir30tif        = './'+outf+'/dir30.tif'
    dir30tau        = './'+outf+'/dir30tau.tif'
    dir30taud4      = './'+outf+'/dir30taud4.tif'
    dir30tau_mask   = './'+outf+'/dir30tau_mask.tif'
    dir30tau_maskd4 = './'+outf+'/dir30tau_maskd4.tif'
    acc3tif         = './'+outf+'/acc3.tif'
    acc30tif        = './'+outf+'/acc30.tif'
    net3tif         = './'+outf+'/net3.tif'
    net30tif        = './'+outf+'/net30.tif'
    net3tifd4       = './'+outf+'/net3d4.tif'
    net30tifd4      = './'+outf+'/net30d4.tif'
    strn_ord3d8     = './'+outf+'/strn_ord3d8.tif'
    strn_tree3d8    = './'+outf+'/strn_tree3d8.txt'
    strn_coord3d8   = './'+outf+'/strn_coord3d8.txt'
    stren_net3d8    = './'+outf+'/stren_net3d8.out'
    stren_w3d8      = './'+outf+'/stren_w3d8.tif'
    strn_ord3d4     = './'+outf+'/strn_ord3d4.tif'
    strn_tree3d4    = './'+outf+'/strn_tree3d4.txt'
    strn_coord3d4   = './'+outf+'/strn_coord3d4.txt'
    stren_net3d4    = './'+outf+'/stren_net3d4.out'
    stren_w3d4      = './'+outf+'/stren_w3d4.tif'
    strn_ord30d8    = './'+outf+'/strn_ord30d8.tif'
    strn_tree30d8   = './'+outf+'/strn_tree30d8.txt'
    strn_coord30d8  = './'+outf+'/strn_coord30d8.txt'
    stren_net30d8   = './'+outf+'/stren_net30d8.out'
    stren_w30d8     = './'+outf+'/stren_w30d8.tif'
    strn_ord30d4    = './'+outf+'/strn_ord30d4.tif'
    strn_tree30d4   = './'+outf+'/strn_tree30d4.txt'
    strn_coord30d4  = './'+outf+'/strn_coord30d4.txt'
    stren_net30d4   = './'+outf+'/stren_net30d4.out'
    stren_w30d4     = './'+outf+'/stren_w30d4.tif'
    out3shp         = './'+outf+'/out3.shp'
    out3shpd4       = './'+outf+'/out3d4.shp'
    out30shp        = './'+outf+'/out30.shp'
    out30shpd4      = './'+outf+'/out30d4.shp'
    cat3tif         = './'+outf+'/basins3.tif'
    cat30tif        = './'+outf+'/basins30.tif'
    cat30tifd4      = './'+outf+'/basins30d4.tif'
    are3tif         = './'+outf+'/area3.tif'
    are30tif        = './'+outf+'/area30.tif'

    # create list of tiles in the region
    # write_list_files(dem3path,'.tif',dem3list)
    # write_list_files(wth3path,'.tif',wth3list)

    # merge tiles in region, create .vrt files, merge tiles
    # subprocess.call(["gdalbuildvrt","-input_file_list",dem3list,dem3vrt])
    # subprocess.call(["gdalbuildvrt","-input_file_list",wth3list,wth3vrt])

    # clip for region
    # subprocess.call(["gdalwarp","-ot","Float32","-te",str(xmin),str(ymin),str(xmax),str(ymax),"-overwrite","-dstnodata","-9999","-co","BIGTIFF=YES",dem3vrt,dem3tif])
    # subprocess.call(["gdalwarp","-ot","Float32","-te",str(xmin),str(ymin),str(xmax),str(ymax),"-overwrite","-dstnodata","-9999","-co","BIGTIFF=YES",wth3vrt,wth3tif])

    if res == 3:

        # listing files in a text file
        write_list_files(dir3path,'.tif',dir3list)

        # building a .vrt file globally
        subprocess.call(["gdalbuildvrt","-input_file_list",dir3list,dir3vrt])

        # clipping interested area output is a .tif file
        subprocess.call(["gdalwarp","-te",str(xmin),str(ymin),str(xmax),str(ymax),"-overwrite","-co","BIGTIFF=YES",dir3vrt,dir3tif])

        # change directions in dir3 file, TAUDEM accepts only direction ranging from 0-8,
        # Yamazaki's DIR follows 1: east, 2: southeast, 4: south, 8: southwest, 16: west, 32: northwest, 64: north.
        # 128: northeast 0: river mouth, -1: inland depression, -9: undefined (ocean) similar to HydroSHEDS
        # same for dir30 file
        print "converting directions into TAUDEM directions..."
        directions_tau(dir3tif,dir3tau)

        print "calculating area in extent..."
        calculate_area(dir3tau,are3tif)

        # calculating accumulation using TauDEM
        print "calculating flow accumulation via TAUDEM..."
        subprocess.call(["mpiexec","-n","20","aread8","-p",dir3tau,"-ad8",acc3tif,"-nc"])

        # thresholding accumulation to get river network
        print "thresholding accumulation to get river network..."
        rasterthreshold(acc3tif,thresh,'Int16',net3tif)

        # masking directions based on river network
        print "masking directions based on river network..."
        rastermask(dir3tau,net3tif,"Int16",dir3tau_mask)

        print "writing outlets and inland depressions in shapefile..."
        write_outlets(out3shp,dir3tau_mask)

        print "writing basins file..."
        subprocess.call(["gagewatershed","-p",dir3tau,"-gw",cat3tif,"-o",out3shp])

        # running streamnet function fro TauDEM to get "tree" and "coord" files
        if streamnet == 'yes':
            subprocess.call(["mpiexec","-n","10","streamnet","-fel",net3tif,"-p",dir3tau,"-ad8",acc3tif,"-src",net3tif,"-ord",strn_ord3d8,"-tree",strn_tree3d8,"-coord",strn_coord3d8,"-net",stren_net3d8,"-w",stren_w3d8,"-o",out3shp])

        # covnert d8 into d4 directions, produce a d4 river network
        print "creating D4 river network..."
        d82d4(dir3tau_mask,dir3tau_maskd4,net3tifd4)

        print "writing D4 outlets and inland depression in shapefile"
        write_outlets(out3shpd4,dir3tau_maskd4)

        # runnning streamnet to get "tree" and "coord" files for the d4 river network
        if streamnet == 'yes':
            subprocess.call(["mpiexec","-n","10","streamnet","-fel",net3tifd4,"-p",dir3tau_maskd4,"-ad8",acc3tif,"-src",net3tifd4,"-ord",strn_ord3d4,"-tree",strn_tree3d4,"-coord",strn_coord3d4,"-net",stren_net3d4,"-w",stren_w3d4,"-o",out3shpd4])

    # same procedure for 30s resolution
    elif res == 30:

        subprocess.call(["gdalwarp","-te",str(xmin),str(ymin),str(xmax),str(ymax),"-overwrite",dir30merged,dir30tif])

        print "converting directions into TAUDEM directions..."
        directions_tau(dir30tif,dir30tau)

        print "calculating area in extent..."
        calculate_area(dir30tau,are30tif)

        print "calculating flow accumulation via TAUDEM..."
        subprocess.call(["mpiexec","-n","20","aread8","-p",dir30tau,"-ad8",acc30tif,"-nc"])

        print "thresholding accumulation to get river network..."
        rasterthreshold(acc30tif,thresh,'Int16',net30tif)

        print "masking directions based on river network..."
        rastermask(dir30tau,net30tif,"Int16",dir30tau_mask)

        print "writing outlets and inland depressions in shapefile..."
        write_outlets(out30shp,dir30tau_mask)

        print "writing basins file..."
        subprocess.call(["gagewatershed","-p",dir30tau,"-gw",cat30tif,"-o",out30shp])

        if streamnet == 'yes':
            subprocess.call(["mpiexec","-n","10","streamnet","-fel",net30tif,"-p",dir30tau,"-ad8",acc30tif,"-src",net30tif,"-ord",strn_ord30d8,"-tree",strn_tree30d8,"-coord",strn_coord30d8,"-net",stren_net30d8,"-w",stren_w30d8,"-o",out30shp])

        print "creating D4 river network..."
        d82d4(dir30tau_mask,dir30tau_maskd4,net30tifd4)

        print "writing D4 outlets and inland depression in shapefile..."
        write_outlets(out30shpd4,dir30tau_maskd4)

        print "create flow directions map D4..."
        create_dir_d4(dir30taud4,dir30tau,dir30tau_maskd4)

        print "writing basins file D4..."
        subprocess.call(["gagewatershed","-p",dir30taud4,"-gw",cat30tifd4,"-o",out30shpd4])

        if streamnet == 'yes':
            subprocess.call(["mpiexec","-n","10","streamnet","-fel",net30tifd4,"-p",dir30tau_maskd4,"-ad8",acc30tif,"-src",net30tifd4,"-ord",strn_ord30d4,"-tree",strn_tree30d4,"-coord",strn_coord30d4,"-net",stren_net30d4,"-w",stren_w30d4,"-o",out30shpd4])


def directions_tau(inputrast,outputrast):

    """
    Function to use in Shell to change convetion from a DIR file
    HydroSHEDS uses ESRI convention 128,64,32,.. this script
    changes these numbers to TauDEM convention, 1,2,3...
    """

    nodata  = -32768
    data    = get_data(inputrast)
    datageo = get_geo(inputrast)
    datatau = cy_directions_tau(np.int16(data),np.int16(nodata))
    gdalutils.write_raster(np.float64(datatau),outputrast,datageo,"Int16",nodata)


def rasterthreshold(file,thres,fmt,outp):

    """
    Output a raster based on a threshold (larger-equal-than)
    """

    nodata   = -1
    filedata = get_data(file)
    filegeo  = get_geo(file)
    data     = cy_rasterthreshold(np.float64(filedata), np.float64(thres), np.float64(nodata))
    gdalutils.write_raster(np.float64(data),outp,filegeo,fmt,nodata)


def rastermask(file,mask,fmt,outp):

    """
    Mask input array following a defined mask (1,0)
    """

    nodata = -32768
    filedata = get_data(file)
    maskdata = get_data(mask)
    filegeo  = get_geo(file)
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
    dirdata  = get_data(filedir)
    dirgeo   = get_geo(filedir)
    data,net = cy_d82d4(np.int16(dirdata), np.int16(nodata))
    write_raster(np.int16(data),filediro,dirgeo,"Int16",nodata)
    write_raster(np.int16(net),fileneto,dirgeo,"Int16",nodata)


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


def find_neighbours(dat,row,col):

    nei = []
    try:
        nei.append(dat[row,col+1])
    except IndexError:
        nei.append(np.nan)
    try:
        nei.append(dat[row-1,col+1])
    except IndexError:
        nei.append(np.nan)
    try:
        nei.append(dat[row-1,col])
    except IndexError:
        nei.append(np.nan)
    try:
        nei.append(dat[row-1,col-1])
    except IndexError:
        nei.append(np.nan)
    try:
        nei.append(dat[row,col-1])
    except IndexError:
        nei.append(np.nan)
    try:
        nei.append(dat[row+1,col-1])
    except IndexError:
        nei.append(np.nan)
    try:
        nei.append(dat[row+1,col])
    except IndexError:
        nei.append(np.nan)
    try:
        nei.append(dat[row+1,col+1])
    except IndexError:
        nei.append(np.nan)

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
    gdalutils.write_raster(dat,output,geo,"Float32",-9999)

if __name__ == '__main__':
    prepdata(sys.argv[1:])
