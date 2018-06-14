#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: j.sosa@bristol.ac.uk / sosa.jeison@gmail.com

import os
import sys
import getopt
import configparser
import numpy as np
import pandas as pd
import gdalutils
from lfptools import misc_utils

def split(argv):
    
    """
    Info:
    -----
    Split area per basins

    Check if basin area is larger than 100 Km2, if yes, it continues

    Check if number of pixels in river network is larger than 35, if yes, continues

    Clip tree and coord files

    To Clip .tif files, it checks outlet direction and increases 0.1 degrees from 
    boundaries but no in the outlet boundary since it will obstruct the outflow.

    There is a “connections” function which finds connections between links. 
    Basically, creates a NNN_rec.csv file containing all river coordinates, three 
    classifications: 1) classification per “link” (given by TAUDEM), 2) classification 
    per reach (coordinates belonging to longer link) and 3) classification per 
    downstream link (given by TAUDEM), other features were added like strahler number 
    (given by TAUDEM) and distance to outlet (given by TAUDEM).


    Dependencies:
    -------------
    TAUDEM : streamnet, gagewatershed
    GDAL : gdalwarp
    Pandas, Numpy, GDAL Python API, gdalutils, prepdata_utils


    Inputs:
    -------
    basnum : Basin number or the keyword “all” to process
    cattif : A GDAL raster file containing basin numbers (one number per basin)
    demtif : DEM in GDAL raster
    acctif : Accumulation GDAL raster
    nettif : River network GDAL raster
    wthtif : River width GDAL raster
    dirtif : Flow directions GDAL raster
    aretif : Area GDAL raster
    otltif : Outlets GDAL raster
    tretxt : Tree file from TAUDEM
    cootxt : Coord file from TAUDEM
    outdir : Out path


    Outputs (If running at 30s):
    ----------------------------
    NNN_acc.tif
    NNN_coo.csv
    NNN_dem.tif
    NNN_dir.tif
    NNN_net.tif
    NNN_rec.csv
    NNN_tre.csv
    NNN_wth.tif
    """

    opts, args = getopt.getopt(argv,"i:")
    for o, a in opts:
        if o == "-i": inifile  = a
    config = configparser.SafeConfigParser()
    config.read(inifile)

    basnum = str(config.get('split','basnum'))
    cattif = str(config.get('split','cattif'))
    demtif = str(config.get('split','demtif'))
    acctif = str(config.get('split','acctif'))
    nettif = str(config.get('split','nettif'))
    wthtif = str(config.get('split','wthtif'))
    dirtif = str(config.get('split','dirtif'))
    otltif = str(config.get('split','otltif'))
    aretif = str(config.get('split','aretif'))
    tretxt = str(config.get('split','tretxt'))
    cootxt = str(config.get('split','cootxt'))
    outdir = str(config.get('split','outdir'))

    print("    running split.py...")

    # Loading data
    catarr = gdalutils.get_data(cattif)

    # Clip input maps per catchment
    if basnum == "all":
        # Loop over all catchment numbers
        for nc in np.unique(catarr[catarr>0]): # Catchments should be numbered and > 0
            basinsplit(nc,outdir,cattif,demtif,acctif,nettif,wthtif,dirtif,aretif,otltif,tretxt,cootxt)
    else:
        # Process a single catchment
        b = basnum.split(',')
        for nc in np.arange(int(b[0]),int(b[1])):
            basinsplit(nc,outdir,cattif,demtif,acctif,nettif,wthtif,dirtif,aretif,otltif,tretxt,cootxt)

def basinsplit(ncatch,outdir,cattif,demtif,acctif,nettif,wthtif,dirtif,aretif,otltif,tretxt,cootxt):

    # Get extend for every catchment and area
    catarr  = gdalutils.get_data(cattif)

    try:
        dat = catarr==ncatch
    except:
        sys.exit('ERROR invalid basin number')

    catgeo  = gdalutils.get_geo(cattif)
    area    = gdalutils.get_data(aretif)
    outlet  = gdalutils.get_data(otltif)
    direc   = gdalutils.get_data(dirtif)   
    row,col = np.where(dat)
    _sum    = np.sum(dat*area)

    if _sum >= 100: # be sure basin is larger than 100 Km2

        xmin = catgeo[8][min(col)]
        xmax = catgeo[8][max(col)]
        ymin = catgeo[9][max(row)]
        ymax = catgeo[9][min(row)]

        # Clip input rasters
        netarr_tmp,netgeo_tmp = gdalutils.clip_raster(nettif,xmin,ymin,xmax,ymax)
        catarr_tmp,catgeo_tmp = gdalutils.clip_raster(cattif,xmin,ymin,xmax,ymax)

        # Mask only the catchment and fill with zeros
        netarr_tmp = np.where(catarr_tmp==ncatch,netarr_tmp,0)
        
        if netarr_tmp.sum() >= 35:  #be sure river network is long enough

            # Clipping tree and coord files based on nettif > 0, coordinates
            tree  = misc_utils.read_tree_taudem(tretxt)
            coor  = misc_utils.read_coord_taudem(cootxt)
            iy,ix = np.where(netarr_tmp > 0)
            Xrav  = netgeo_tmp[8][ix]
            Yrav  = netgeo_tmp[9][iy]

            # Clipping coord file (it may be improved, calculation takes some time)
            lfp_coor = pd.DataFrame()
            for i in range(len(Xrav)):
                dis,ind = misc_utils.near_euc(coor['lon'].values,coor['lat'].values,(Xrav[i],Yrav[i]))
                if dis <= 0.01:
                    lfp_coor = lfp_coor.append(coor.loc[ind,:])
            lfp_coor = lfp_coor[['lon','lat','distance','elev','contr_area']]
            lfp_coor.index.name = 'index'
            lfp_coor.sort_index(inplace=True)
            lfp_coor.drop_duplicates(inplace=True) # Remove duplicates just in case

            # Clipping tree file
            lfp_tree = pd.DataFrame()
            for i in tree.index:
                sta  = tree.loc[i,'start_pnt']
                end  = tree.loc[i,'end_pnt']
                lon1 = coor.loc[sta,'lon']
                lat1 = coor.loc[sta,'lat']
                lon2 = coor.loc[end,'lon']
                lat2 = coor.loc[end,'lat']
                dis1,ind1 = misc_utils.near_euc(lfp_coor['lon'].values,lfp_coor['lat'].values,(lon1,lat1))
                dis2,ind2 = misc_utils.near_euc(lfp_coor['lon'].values,lfp_coor['lat'].values,(lon2,lat2))
                if (dis1 <= 0.012) & (dis2 <= 0.012): # default value 0.01 wasn't able to find link number 3504, this value was increased to 0.012 to find missing link
                    lfp_tree = lfp_tree.append(tree.loc[i,:])
            lfp_tree = lfp_tree[['link_no','start_pnt','end_pnt','frst_ds','frst_us','scnd_us','strahler','mon_pnt','shreve']]
            lfp_tree.index.name = 'index'

            # Creating folder per basin
            ncatchstr = "%03d" % ncatch
            folder = outdir + "/" + ncatchstr
            create_out_folder(folder)

            # Writing clipped coord and tree files
            fnametre = folder + "/" + ncatchstr +"_tre.csv"
            fnamecoo = folder + "/" + ncatchstr +"_coo.csv"
            lfp_coor.to_csv(fnamecoo)
            lfp_tree.to_csv(fnametre,float_format='%i')

            # Creating rec dataframe
            rec = connections(fnametre,fnamecoo)

            #  Writing XXX_rec.csv file
            fnamerec = folder + "/" + ncatchstr +"_rec.csv"
            rec.to_csv(fnamerec)

            # Get extent from rec dataframe
            xmin = rec['lon'].min()
            xmax = rec['lon'].max()
            ymin = rec['lat'].min()
            ymax = rec['lat'].max()

            # Get fixed extent
            # _dir    = getdir(rec,dirtif)
            # _dirlet = getdirletter(_dir)
            # xmin,ymin,xmax,ymax = get_extent_outlet(_dirlet,0.1,xmin,ymin,xmax,ymax)

            # Clipping rasters
            demarrcli,demgeocli = gdalutils.clip_raster(demtif,xmin,ymin,xmax,ymax)
            accarrcli,accgeocli = gdalutils.clip_raster(acctif,xmin,ymin,xmax,ymax)
            wtharrcli,wthgeocli = gdalutils.clip_raster(wthtif,xmin,ymin,xmax,ymax)
            dirarrcli,dirgeocli = gdalutils.clip_raster(dirtif,xmin,ymin,xmax,ymax)
            netarrcli,netgeocli = gdalutils.clip_raster(nettif,xmin,ymin,xmax,ymax)
            catarrcli,catgeocli = gdalutils.clip_raster(cattif,xmin,ymin,xmax,ymax)

            # Mask only the catchment and fill with zeros
            netarrcli = np.where(catarrcli==ncatch,netarrcli,0)
            dirarrcli = np.where(catarrcli==ncatch,dirarrcli,0)

            # Creating output names
            fnamedem = folder + "/" + ncatchstr +"_dem.tif"
            fnameacc = folder + "/" + ncatchstr +"_acc.tif"
            fnamenet = folder + "/" + ncatchstr +"_net.tif"
            fnamewth = folder + "/" + ncatchstr +"_wth.tif"
            fnamedir = folder + "/" + ncatchstr +"_dir.tif"

            # Writing clipped arrays
            nodata = -9999
            gdalutils.write_raster(demarrcli,fnamedem,demgeocli,"Float32",nodata)
            gdalutils.write_raster(accarrcli,fnameacc,accgeocli,"Float32",nodata)
            gdalutils.write_raster(netarrcli,fnamenet,netgeocli,"Float32",nodata)
            gdalutils.write_raster(wtharrcli,fnamewth,wthgeocli,"Float32",nodata)
            gdalutils.write_raster(dirarrcli,fnamedir,dirgeocli,"Float32",nodata)

        else:
            print("NOT PROCESSED: Number of pixels in river lower than 35 : " + str(netarr_tmp.sum()) + " pixels in basin number " + str(ncatch))
    else:
        print("NOT PROCESSED: Basin area lower than 100 Km2 : " + str(_sum) + " KM**2 in basin number " + str(ncatch))

def connections(treef,coorf):
    
    """
    Finds connections between links and sort them
    First finds the connections for all links, then sort them
    starting from the link with more downstream connections, later
    write separate files including coordinates and links
    """
    
    def find_links(link):
        mylinks =[]
        mylinks.append(link)
        while True:
            linkds = tree.loc[link,'frst_ds']
            if linkds == -1: break
            else:
                mylinks.append(linkds)
                link = linkds
        return mylinks
    
    tree = misc_utils.read_tree(treef)
    coor = misc_utils.read_coord(coorf)
    tree.set_index('link_no',inplace=True)
    
    # Finding the number of downstream links for every link
    lnks = []
    size = []
    for i in tree.index:
        res = find_links(i)
        lnks.append(res)
        size.append(len(res))

    # Create columns with number of downstream links, index and flag with zeros
    tree['links'] = size
    tree['index'] = range(tree.index.size)
    tree['link_flag'] = 0

    # Sorting in a descending way, links with more downstream links go first
    tree.sort_values(by='links',ascending=False,inplace=True)

    # Go over links, check links from the lnks list
    c = 0
    df_rec = pd.DataFrame()
    for i in tree['index']:
        df_cor = pd.DataFrame()
        for j in lnks[i]:
            if tree.loc[j,'link_flag'] == 0:
                tree.loc[j,'link_flag'] = 1
                # Retrieve elevation
                start = tree.loc[j,'start_pnt']
                end = tree.loc[j,'end_pnt']
                df = coor.loc[start:end,'lon':'distance']
                df['link'] = int(j)
                df_cor = pd.concat([df_cor,df])
        if df_cor.empty is False:
            c+=1
            df_cor['reach'] = int(c)
            df_rec = pd.concat([df_rec,df_cor])

    # Retrieving Strahler number
    dslk = []
    stra = []
    for i in df_rec.index:
        link = df_rec.loc[i,'link']
        stra_val = tree.loc[link,'strahler']
        dslk_val = tree.loc[link,'frst_ds']
        stra.append(stra_val)
        dslk.append(dslk_val)
    df_rec['strahler'] = stra
    df_rec['dslink'] = dslk

    return df_rec

def create_out_folder(folder):

    try:
        os.makedirs(folder)
    except OSError:
        if not os.path.isdir(folder):
            raise

def getdirletter(dirval):

    if dirval == 1:
        dirlet = "E"
    elif dirval == 3:
        dirlet = "N"
    elif dirval == 5:
        dirlet = "W"
    elif dirval == 7:
        dirlet = "S"
    else:
        sys.exit('ERROR: Wrong direction found')
    return dirlet

def get_extent_outlet(dirlet,thresh,_xmin,_ymin,_xmax,_ymax):

    if dirlet == "E":

        xmin = _xmin - thresh
        ymin = _ymin - thresh
        xmax = _xmax
        ymax = _ymax + thresh

    elif dirlet == "W":

        xmin = _xmin
        ymin = _ymin - thresh
        xmax = _xmax + thresh
        ymax = _ymax + thresh

    elif dirlet == "N":

        xmin = _xmin - thresh
        ymin = _ymin - thresh
        xmax = _xmax + thresh
        ymax = _ymax
        
    elif dirlet == "S":

        xmin = _xmin - thresh
        ymin = _ymin
        xmax = _xmax + thresh
        ymax = _ymax + thresh   

    return (xmin,ymin,xmax,ymax)

def getdir(rec,dirtif):
    
    dat   = gdalutils.get_data(dirtif)
    geo   = gdalutils.get_geo(dirtif)
    dirdf = gdalutils.array_to_pandas(dat,geo,0,'gt')
    recdf = gdalutils.assign_val(df2=rec.reset_index(),df2_x='lon',df2_y='lat',df1=dirdf,df1_x='x',df1_y='y',label='z',copy=True)

    # Direction of outlet is given by the maximum repetitions of directions in the last 10 points
    _dir  = recdf.sort_values(by='distance').iloc[0:10].groupby('z')['z'].count().idxmax()

    return _dir

if __name__ == '__main__':
    split(sys.argv[1:])
