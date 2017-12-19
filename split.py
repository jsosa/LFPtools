#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# date: 26/jun/2017
# mail: j.sosa@bristol.ac.uk / sosa.jeison@gmail.com

import os,sys,getopt
import ConfigParser
from gdal_utils import *
from scipy.spatial.distance import cdist
from scipy.ndimage import distance_transform_edt

import pdb
# pdb.set_trace()

def split(argv):

    opts, args = getopt.getopt(argv,"i:")
    for o, a in opts:
        if o == "-i": inifile  = a
    config = ConfigParser.SafeConfigParser()
    config.read(inifile)

    basinnum  = str(config.get('split','basinnum'))
    dem3fil   = str(config.get('split','dem3fil'))
    dir3fil   = str(config.get('split','dir3fil'))
    acc3fil   = str(config.get('split','acc3fil'))
    net3fil   = str(config.get('split','net3fil'))
    net30fil  = str(config.get('split','net30fil'))
    dir30fil  = str(config.get('split','dir30fil'))
    tre30fil  = str(config.get('split','tre30fil'))
    coo30fil  = str(config.get('split','coo30fil'))
    cat30fil  = str(config.get('split','cat30fil'))
    outdirfil = str(config.get('split','outdirfil'))

    # loading data
    cat30geo = get_gdal_geo(cat30fil)
    cat30dat = get_gdal_data(cat30fil)

    ### CLIP INPUT MAPS PER CATCHMENT
    
    if basinnum == "all":

        # loop over all catchment numbers
        for nc in np.unique(cat30dat[cat30dat>0]): # catchments should be numbered >0
            print "split.py - " + str(np.max(cat30dat)-nc)
            basinsplit(nc,outdirfil,cat30geo,cat30dat,dem3fil,acc3fil,net3fil,net30fil,dir3fil,dir30fil,cat30fil,coo30fil,tre30fil)
    else:
        # process a single catchment
        print "split.py - 1"
        basinsplit(int(basinnum),outdirfil,cat30geo,cat30dat,dem3fil,acc3fil,net3fil,net30fil,dir3fil,dir30fil,cat30fil,coo30fil,tre30fil)

def basinsplit(ncatch,outdirfil,cat30geo,cat30dat,dem3fil,acc3fil,net3fil,net30fil,dir3fil,dir30fil,cat30fil,coo30fil,tre30fil):

    ncatchstr = "%03d" % ncatch

    # create basin folder
    folder = outdirfil + "/" + ncatchstr
    try:
        os.makedirs(folder)
    except OSError:
        if not os.path.isdir(folder):
            raise

    row,col = np.where(cat30dat==ncatch) # get two vectors specifying index where clausule is valid

    # This threshold produce an extra pixel after the outlet of the river!
    # xmin = cat30geo[8][min(col)] - 0.01
    # xmax = cat30geo[8][max(col)] + 0.01
    # ymin = cat30geo[9][max(row)] - 0.01
    # ymax = cat30geo[9][min(row)] + 0.01

    xmin = cat30geo[8][min(col)]
    xmax = cat30geo[8][max(col)]
    ymin = cat30geo[9][max(row)]
    ymax = cat30geo[9][min(row)]

    dem3datcli,dem3geocli   = clip_raster(dem3fil,xmin,ymin,xmax,ymax)
    acc3datcli,acc3geocli   = clip_raster(acc3fil,xmin,ymin,xmax,ymax)
    dir3datcli,dir3geocli   = clip_raster(dir3fil,xmin,ymin,xmax,ymax)
    net3datcli,net3geocli   = clip_raster(net3fil,xmin,ymin,xmax,ymax)
    net30datcli,net30geocli = clip_raster(net30fil,xmin,ymin,xmax,ymax)
    dir30datcli,dir30geocli = clip_raster(dir30fil,xmin,ymin,xmax,ymax)
    cat30datcli,cat30geocli = clip_raster(cat30fil,xmin,ymin,xmax,ymax)

    # mask only the catchment

    nodata = -9999
    
    net30datcli = np.where(cat30datcli==ncatch,net30datcli,0)
    dir30datcli = np.where(cat30datcli==ncatch,dir30datcli,0)

    fnamedem3  = folder + "/" + ncatchstr +"_dem3.tif"
    fnameacc3  = folder + "/" + ncatchstr +"_acc3.tif"
    fnamedir3  = folder + "/" + ncatchstr +"_dir3.tif"
    fnamenet3  = folder + "/" + ncatchstr +"_net3.tif"
    fnamenet30 = folder + "/" + ncatchstr +"_net30.tif"
    fnamedir30 = folder + "/" + ncatchstr +"_dir30.tif"

    writeRaster(dem3datcli,fnamedem3,dem3geocli,"Float32",nodata)
    writeRaster(acc3datcli,fnameacc3,acc3geocli,"Float32",nodata)
    writeRaster(dir3datcli,fnamedir3,dir3geocli,"Float32",nodata)
    writeRaster(net3datcli,fnamenet3,net3geocli,"Float32",nodata)
    writeRaster(net30datcli,fnamenet30,net30geocli,"Float32",nodata)
    writeRaster(dir30datcli,fnamedir30,cat30geocli,"Float32",nodata)

    ### WRITE TREE.TXT FILE PER CATCHMENT ###

    coord = np.genfromtxt(coo30fil, delimiter="\t")
    coord = np.delete(coord,0,1) # remove first column, is an empty column
    tree  = np.genfromtxt(tre30fil, delimiter="\t")
    tree  = np.delete(tree,0,1) # remove first column, is an empty column

    file = open(outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_tre.txt","w")

    nlinks = 0
    for i in range(tree.shape[0]):
        
        start = int(tree[i,1])
        end   = int(tree[i,2])
        lon1  = coord[start,0]
        lat1  = coord[start,1]
        lon2  = coord[end,0]
        lat2  = coord[end,1]

        indy,indx = np.where(net30datcli>0)
        degx,degy = net30geocli[8][indx],net30geocli[9][indy]

        if near(degx,degy,np.array([[lat1,lon1]])) < 0.01 and near(degx,degy,np.array([[lat2,lon2]])) < 0.01:
            nlinks = nlinks + 1
            file.write("%d"%tree[i,0]+" "+"%d"%tree[i,1]+" "+"%d"%tree[i,2]+" "+"%d"%tree[i,3]+" "+"%d"%tree[i,4]+"\n")

    file.close()

    ### WRITE COORD.TXT FILE PER CATCHMENT ###

    # Check if tree links are inside the basin, this validation is required since we use the
    # Euclidean distance of "0.1" to determine if a link is inside (from previous strep)
    tree_cli = np.genfromtxt(outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_tre.txt", delimiter=" ")
    file = open(outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_coo.txt","w")

    file.write("lon"+" "+"lat"+"\n")
    for itree in range(nlinks):

        if nlinks == 1:
            start = int(tree_cli[1])
            end   = int(tree_cli[2]) + 1 # inclusive index
        else:
            start = int(tree_cli[itree,1])
            end   = int(tree_cli[itree,2]) + 1 # inclusive index

        lons  = coord[start:end,0]
        lats  = coord[start:end,1]
        
        for ilons in range(lons.shape[0]):
            file.write("%.8f"%lons[ilons]+" "+"%.8f"%lats[ilons]+"\n")
    file.close()

    ### WRITE INI FILES PER CATCHMENT ###

    config = ConfigParser.RawConfigParser()
    
    # getinflows.py
    config.add_section('getinflows')
    config.set('getinflows', 'thresh', '4')
    config.set('getinflows', 'output',  outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_inf")
    config.set('getinflows', 'method', 'haversine')
    config.set('getinflows', 'netf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_net30.tif")
    config.set('getinflows', 'treef',   outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_tre.txt")
    config.set('getinflows', 'coordf', '/Users/js16170/SOSA001/university_of_bristol/js16170/projects/euflood/dat/eu_coordd430s.txt')
    config.set('getinflows', 'proj',   '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

    # getwidths.py
    config.add_section('getwidths')
    config.set('getwidths', 'thresh', '0.01')
    config.set('getwidths', 'output',  outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_wdt")
    config.set('getwidths', 'method', 'near')
    config.set('getwidths', 'netf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_net30.tif")
    config.set('getwidths', 'proj',   '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    config.set('getwidths', 'fwidth', '/Users/js16170/SOSA001/university_of_bristol/js16170/projects/euflood/dat/eu_width_andreadis.tif')

    # getbankelevs.py
    config.add_section('getbankelevs')
    config.set('getbankelevs', 'output',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_bnk")
    config.set('getbankelevs', 'netf',      outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_net30.tif")
    config.set('getbankelevs', 'hrdemf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_dem3.tif")
    config.set('getbankelevs', 'hrrivf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_net3.tif")
    config.set('getbankelevs', 'proj',     '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    config.set('getbankelevs', 'method',   'mean')
    config.set('getbankelevs', 'outlier',  'yes')
    config.set('getbankelevs', 'hrnodata', '-9999')
    config.set('getbankelevs', 'thresh',   '0.00416')

    # fixelevs.py
    config.add_section('fixelevs')
    config.set('fixelevs', 'source',  outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_bnk.shp")
    config.set('fixelevs', 'output',  outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_bnkfix")
    config.set('fixelevs', 'netf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_net30.tif")
    config.set('fixelevs', 'coordf', '/Users/js16170/SOSA001/university_of_bristol/js16170/projects/euflood/dat/eu_coordd430s.txt')
    config.set('fixelevs', 'treef',   outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_tre.txt")
    config.set('fixelevs', 'proj',   '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    config.set('fixelevs', 'method', 'yamazaki')

    # getslopes.py
    config.add_section('getslopes')
    config.set('getslopes', 'source',  outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_bnkfix.shp")
    config.set('getslopes', 'output',  outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_slopes")
    config.set('getslopes', 'netf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_net30.tif")
    config.set('getslopes', 'coordf', '/Users/js16170/SOSA001/university_of_bristol/js16170/projects/euflood/dat/eu_coordd430s.txt')
    config.set('getslopes', 'treef',   outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_tre.txt")
    config.set('getslopes', 'proj',   '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    config.set('getslopes', 'method', 'linear_regr_step')
    config.set('getslopes', 'step',   '5')

    # rasterresample.py
    config.add_section('rasterresample')
    config.set('rasterresample', 'method',   'mean')
    config.set('rasterresample', 'demf',      outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_dem3.tif")
    config.set('rasterresample', 'netf',      outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_net30.tif")
    config.set('rasterresample', 'output',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_dem30.tif")
    config.set('rasterresample', 'outlier',  'yes')
    config.set('rasterresample', 'hrnodata', '-9999')
    config.set('rasterresample', 'thresh',   '0.00416')
    config.set('rasterresample', 'resx',     '0.0083333333333')
    config.set('rasterresample', 'resy',     '0.0083333333333')
    config.set('rasterresample', 'nproc',    '4') # number of cpus to use

    # getqbankfull.py
    config.add_section('getqbankfull')
    config.set('getqbankfull', 'netf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_net30.tif")
    config.set('getqbankfull', 'output',  outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_qbfull")
    config.set('getqbankfull', 'proj',   '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    config.set('getqbankfull', 'thresh', '5000') # given in meters since JRC data is in meters!
    config.set('getqbankfull', 'aep',    '1.5') # Annual excende probability
    config.set('getqbankfull', 'nproc',  '4') # number of cpus to use

    # getdepths.py
    config.add_section('getdepths')
    config.set('getdepths', 'proj',   '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    config.set('getdepths', 'netf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_net30.tif")
    config.set('getdepths', 'method', 'depth_manning')
    config.set('getdepths', 'output',  outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_dpt")
    config.set('getdepths', 'n',      '0.025')
    config.set('getdepths', 'wdtf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_wdt.shp")
    config.set('getdepths', 'slpf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_slopes.shp")
    config.set('getdepths', 'qbnkf',   outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_qbfull.shp")

    # getbedelevs.py
    config.add_section('getbedelevs')
    config.set('getbedelevs', 'bnkf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_bnkfix.shp")
    config.set('getbedelevs', 'dptf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_dpt.shp")
    config.set('getbedelevs', 'netf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_net30.tif")
    config.set('getbedelevs', 'output',  outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_bed")
    config.set('getbedelevs', 'proj',   '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

    # wirtefiles.py
    config.add_section('writefiles')
    config.set('writefiles', 'parf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+".par")
    config.set('writefiles', 'demf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_dem30.tif")
    config.set('writefiles', 'grpf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_grp.tif")
    config.set('writefiles', 'bnkf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_bnk.tif")
    config.set('writefiles', 'wdtf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_wdt.tif")
    config.set('writefiles', 'inff',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_inf.tif")
    config.set('writefiles', 'bedf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_bed.tif")
    config.set('writefiles', 'fixbnkf', outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_bnkfix.tif")
    config.set('writefiles', 'dembnkf', outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_dembnk.tif")
    config.set('writefiles', 'grdcf',   outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_grdc.txt")
    config.set('writefiles', 'nrfaf',   outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_nrfa.txt")
    config.set('writefiles', 'pramf',   outdirfil+"/"+ncatchstr+"/"+ncatchstr+".pram")
    config.set('writefiles', 'netf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_net30.tif")
    config.set('writefiles', 'gaugef',  outdirfil+"/"+ncatchstr+"/"+ncatchstr+".gauge")
    config.set('writefiles', 'stagef',  outdirfil+"/"+ncatchstr+"/"+ncatchstr+".stage")
    config.set('writefiles', 'dirf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_dir30.tif")
    config.set('writefiles', 'evapf',   outdirfil+"/"+ncatchstr+"/"+ncatchstr+".evap")
    config.set('writefiles', 'bdyf',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+".bdy")
    config.set('writefiles', 'bcif',    outdirfil+"/"+ncatchstr+"/"+ncatchstr+".bci")
    config.set('writefiles', 'date1',   '1990-01-01')
    config.set('writefiles', 'date2',   '1990-01-10')
    
    with open(outdirfil+"/"+ncatchstr+"/config.cfg", 'wb') as configfile:
        config.write(configfile)

    ### WRITE main.sh FILE ###

    toolbox = '/Users/js16170/SOSA001/university_of_bristol/js16170/projects/euflood/src/'

    file = open(outdirfil+"/"+ncatchstr+"/"+ncatchstr+"_main.sh","w")
    file.write('python'+' '+toolbox+"getinflows.py"+' '+'-i'+' '+'config.cfg'+'\n')
    file.write('python'+' '+toolbox+"getwidths.py"+' '+'-i'+' '+'config.cfg'+'\n')
    file.write('python'+' '+toolbox+"getbankelevs.py"+' '+'-i'+' '+'config.cfg'+'\n')
    file.write('python'+' '+toolbox+"fixelevs.py"+' '+'-i'+' '+'config.cfg'+'\n')
    file.write('python'+' '+toolbox+"getslopes.py"+' '+'-i'+' '+'config.cfg'+'\n')
    file.write('python'+' '+toolbox+"getqbankfull.py"+' '+'-i'+' '+'config.cfg'+'\n')
    file.write('python'+' '+toolbox+"getdepths.py"+' '+'-i'+' '+'config.cfg'+'\n')
    file.write('python'+' '+toolbox+"getbedelevs.py"+' '+'-i'+' '+'config.cfg'+'\n')
    file.write('python'+' '+toolbox+"rasterresample.py"+' '+'-i'+' '+'config.cfg'+'\n')
    file.write('python'+' '+toolbox+"writefiles.py"+' '+'-i'+' '+'config.cfg'+'\n')
    file.close()

def near(ddsx,ddsy,XA):

    XB  = np.vstack((ddsy,ddsx)).T
    dis = cdist(XA, XB, metric='euclidean').min()

    return dis

if __name__ == '__main__':
    split(sys.argv[1:])