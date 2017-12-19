#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# date: 22/nov/2017
# mail: j.sosa@bristol.ac.uk / sosa.jeison@gmail.com

import os
import sys
import pdb # pdb.set_trace()
import getopt
import subprocess
import numpy as np
from osgeo import osr
from osgeo import gdal
from prepdata_utils import cy_d82d4
from prepdata_utils import cy_rastermask
from prepdata_utils import cy_directions_tau
from prepdata_utils import cy_rasterthreshold

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

    wth3path = '/scratch/hydro3/js16170/data_processed/yamazaki/wth3/'

    # list of files generated
    dem3list        = './'+outf+'/dem3.txt'
    dem3vrt         = './'+outf+'/dem3.vrt'
    dem3tif         = './'+outf+'/dem3.tif'
    wth3list        = './'+outf+'/wth3.txt'
    wth3vrt         = './'+outf+'/wht3.vrt'
    wth3tif         = './'+outf+'/wht3.tif'
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

        # calculating accumulation using TauDEM
        subprocess.call(["mpiexec","-n","20","aread8","-p",dir3tau,"-ad8",acc3tif,"-nc"])

        # thresholding accumulation to get river network
        print "thresholding accumulation to get river network..."
        rasterthreshold(acc3tif,thresh,'Int16',net3tif)

        # running streamnet function fro TauDEM to get "tree" and "coord" files
        if streamnet == 'yes':
            subprocess.call(["mpiexec","-n","10","streamnet","-fel",net3tif,"-p",dir3tau,"-ad8",acc3tif,"-src",net3tif,"-ord",strn_ord3d8,"-tree",strn_tree3d8,"-coord",strn_coord3d8,"-net",stren_net3d8,"-w",stren_w3d8])

        # masking directions based on river network
        print "masking directions based on river network..."
        rastermask(dir3tau,net3tif,"Int16",dir3tau_mask)

        # covnert d8 into d4 directions, produce a d4 river network
        print "creating D4 river network..."
        d82d4(dir3tau_mask,dir3tau_maskd4,net3tifd4)

        # runnning streamnet to get "tree" and "coord" files for the d4 river network
        if streamnet == 'yes':
            subprocess.call(["mpiexec","-n","10","streamnet","-fel",net3tifd4,"-p",dir3tau_maskd4,"-ad8",acc3tif,"-src",net3tifd4,"-ord",strn_ord3d4,"-tree",strn_tree3d4,"-coord",strn_coord3d4,"-net",stren_net3d4,"-w",stren_w3d4])

    # same procedure for 30s resolution
    elif res == 30:

        write_list_files(dir30path,'.bil',dir30list)

        subprocess.call(["gdalbuildvrt","-input_file_list",dir30list,dir30vrt])

        subprocess.call(["gdalwarp","-te",str(xmin),str(ymin),str(xmax),str(ymax),"-overwrite",dir30vrt,dir30tif])

        print "converting directions into TAUDEM directions..."
        directions_tau(dir30tif,dir30tau)

        subprocess.call(["mpiexec","-n","20","aread8","-p",dir30tau,"-ad8",acc30tif,"-nc"])

        print "thresholding accumulation to get river network..."
        rasterthreshold(acc30tif,thresh,'Int16',net30tif)

        if streamnet == 'yes':
            subprocess.call(["mpiexec","-n","10","streamnet","-fel",net30tif,"-p",dir30tau,"-ad8",acc30tif,"-src",net30tif,"-ord",strn_ord30d8,"-tree",strn_tree30d8,"-coord",strn_coord30d8,"-net",stren_net30d8,"-w",stren_w30d8])

        print "masking directions based on river network..."
        rastermask(dir30tau,net30tif,"Int16",dir30tau_mask)

        print "creating D4 river network..."
        d82d4(dir30tau_mask,dir30tau_maskd4,net30tifd4)

        if streamnet == 'yes':
            subprocess.call(["mpiexec","-n","10","streamnet","-fel",net30tifd4,"-p",dir30tau_maskd4,"-ad8",acc30tif,"-src",net30tifd4,"-ord",strn_ord30d4,"-tree",strn_tree30d4,"-coord",strn_coord30d4,"-net",stren_net30d4,"-w",stren_w30d4])

def get_gdal_data(filename):

    """ Import gdal data in a numpy array """

    ds   = gdal.Open(filename, gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    data = np.float64(band.ReadAsArray())
    ds   = None
    return data

def get_gdal_geo(filename):

    """ Read geo info from raster file """

    ds   = gdal.Open(filename, gdal.GA_ReadOnly)
    gt   = ds.GetGeoTransform()
    nx   = np.int64(ds.RasterXSize)
    ny   = np.int64(ds.RasterYSize)
    resx = np.float64(gt[1])
    resy = np.float64(gt[5])
    xmin = np.float64(gt[0])
    ymin = np.float64(gt[3]) + ny * resy
    xmax = np.float64(gt[0]) + nx * resx
    ymax = np.float64(gt[3])
    x    = np.linspace(xmin+resx/2,xmin+resx/2+resx*(nx-1),nx)
    y    = np.linspace(ymax+resy/2,ymax+resy/2+resy*(ny-1),ny)
    srs  = osr.SpatialReference()
    wkt  = ds.GetProjectionRef()
    srs.ImportFromWkt(wkt)
    ds   = None
    geo  = [xmin,ymin,xmax,ymax,nx,ny,resx,resy,x,y,srs]
    return geo

def writeRaster(myarray,myraster,geo,fmt,nodata):

    """ Write an array to Raster format defaul in GTiff format"""

    # available data types:
    if fmt == "Byte":
        ofmt = gdal.GDT_Byte
    elif fmt == "Float32":
        ofmt = gdal.GDT_Float32
    elif fmt == "Float64":
        ofmt = gdal.GDT_Float64
    elif fmt == "Int16":
        ofmt = gdal.GDT_Int16
    elif fmt == "Int32":
        ofmt = gdal.GDT_Int32

    xmin = geo[0]
    ymin = geo[1]
    xmax = geo[2]
    ymax = geo[3]
    nx   = geo[4]
    ny   = geo[5]
    resx = geo[6]
    resy = geo[7]
    x    = geo[8]
    y    = geo[9]
    srs  = geo[10]

    driver = gdal.GetDriverByName('GTiff')

    outRaster = driver.Create(myraster, nx, ny, 1, ofmt, ['COMPRESS=LZW'])
    outRaster.SetGeoTransform((xmin, resx, 0, ymax, 0, resy))
    outRaster.SetProjection(srs.ExportToWkt())

    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(myarray)
    outband.FlushCache()
    outband.SetNoDataValue(nodata)

def directions_tau(inputrast,outputrast):

    """
    Function to use in Shell to change convetion from a DIR file
    HydroSHEDS uses ESRI convention 128,64,32,.. this script
    changes these numbers to TauDEM convention, 1,2,3...
    """

    nodata  = -32768
    data    = get_gdal_data(inputrast)
    datageo = get_gdal_geo(inputrast)
    datatau = cy_directions_tau(np.int16(data),np.int16(nodata))
    writeRaster(np.float64(datatau),outputrast,datageo,"Int16",nodata)

def rasterthreshold(file,thres,fmt,outp):

    """
    Output a raster based on a threshold (larger-equal-than)
    """

    nodata   = -1
    filedata = get_gdal_data(file)
    filegeo  = get_gdal_geo(file)
    data     = cy_rasterthreshold(np.float64(filedata), np.float64(thres), np.float64(nodata))
    writeRaster(np.float64(data),outp,filegeo,fmt,nodata)

def rastermask(file,mask,fmt,outp):

    """
    Mask input array following a defined mask (1,0)
    """

    nodata = -32768
    filedata = get_gdal_data(file)
    maskdata = get_gdal_data(mask)
    filegeo  = get_gdal_geo(file)
    data     = cy_rastermask(np.float64(filedata),np.int16(maskdata))
    writeRaster(np.float64(data),outp,filegeo,fmt,nodata)

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
    dirdata  = get_gdal_data(filedir)
    dirgeo   = get_gdal_geo(filedir)
    data,net = cy_d82d4(np.int16(dirdata), np.int16(nodata))
    writeRaster(np.int16(data),filediro,dirgeo,"Int16",nodata)
    writeRaster(np.int16(net),fileneto,dirgeo,"Int16",nodata)

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

def outflows_from_dirs(srcdir,foutput):

    mytiles = listdir(srcdir,'.tif')

    w = sf.Writer(sf.POINT)
    w.field('x')
    w.field('y')

    for i in range(len(mytiles)):
        print mytiles[i]
        dat   = get_gdal_data(mytiles[i])
        geo   = get_gdal_geo(mytiles[i])
        iy,ix = np.where(dat==0)
        x     = geo[8][ix]
        y     = geo[9][iy]

        if len(x)>0:
            for j in range(len(x)):
                w.point(x[j],y[j])
                w.record(x[j],y[j])

    w.save("%s.shp" % foutput)

    proj = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    # write .prj file
    prj = open("%s.prj" % foutput, "w")
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj)
    prj.write(srs.ExportToWkt())
    prj.close()

if __name__ == '__main__':
    prepdata(sys.argv[1:])
