#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# date: 26/jun/2017
# mail: j.sosa@bristol.ac.uk / sosa.jeison@gmail.com

import os
import sys
import getopt
import ConfigParser
import numpy as np
import pandas as pd
import gdal_utils
import misc_utils

def split(argv):
    
    """
    Split area by catchments

    Parameters
    ----------

    basnum : basin number
    cattif : catchment file masking every catchment
    demtif : DEM
    acctif : accumulation
    nettif : river network
    wthtif : river width
    dirtif : direction
    tretxt : tree file
    cootxt : coordinates file
    outdir : output directory

    Returns
    -------

    XXX_dem.tif : DEM
    XXX_net.tif : river network mask
    XXX_wth.tif : river width
    XXX_acc.tif : accumulation
    XXX_dir.tif : direction
    XXX_tre.csv : tree file
    XXX_coo.csv : coordinates dile
    XXX_rec.csv : summary file
    """

    opts, args = getopt.getopt(argv,"i:")
    for o, a in opts:
        if o == "-i": inifile  = a
    config = ConfigParser.SafeConfigParser()
    config.read(inifile)

    basnum = str(config.get('split','basnum'))
    cattif = str(config.get('split','cattif'))
    demtif = str(config.get('split','demtif'))
    acctif = str(config.get('split','acctif'))
    nettif = str(config.get('split','nettif'))
    wthtif = str(config.get('split','wthtif'))
    dirtif = str(config.get('split','dirtif'))
    tretxt = str(config.get('split','tretxt'))
    cootxt = str(config.get('split','cootxt'))
    outdir = str(config.get('split','outdir'))

    # Loading data
    catarr = gdal_utils.get_gdal_data(cattif)

    # Clip input maps per catchment
    if basnum == "all":
        # Loop over all catchment numbers
        for nc in np.unique(catarr[catarr>0]): # Catchments should be numbered and > 0
            print "split.py - " + str(np.max(catarr)-nc)
    else:
        # Process a single catchment
        print "split.py - " + basnum
        basinsplit(np.int(basnum),outdir,cattif,demtif,acctif,nettif,wthtif,dirtif,tretxt,cootxt)


def basinsplit(ncatch,outdir,cattif,demtif,acctif,nettif,wthtif,dirtif,tretxt,cootxt):
    
    # Create basin folder
    ncatchstr = "%03d" % ncatch
    folder = outdir + "/" + ncatchstr
    try:
        os.makedirs(folder)
    except OSError:
        if not os.path.isdir(folder):
            raise

    # Get extend for every catchment
    catarr  = gdal_utils.get_gdal_data(cattif)
    catgeo  = gdal_utils.get_gdal_geo(cattif)
    row,col = np.where(catarr==ncatch)
    xmin    = catgeo[8][min(col)]
    xmax    = catgeo[8][max(col)]
    ymin    = catgeo[9][max(row)]
    ymax    = catgeo[9][min(row)]

    # Clip input rasters
    netarr_tmp,netgeo_tmp = gdal_utils.clip_raster(nettif,xmin,ymin,xmax,ymax)
    catarr_tmp,catgeo_tmp = gdal_utils.clip_raster(cattif,xmin,ymin,xmax,ymax)

    # Mask only the catchment and fill with zeros
    netarr_tmp = np.where(catarr_tmp==ncatch,netarr_tmp,0)

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
        sta = tree.loc[i,'start_pnt']
        end = tree.loc[i,'end_pnt']
        lon1 = coor.loc[sta,'lon']
        lat1 = coor.loc[sta,'lat']
        lon2 = coor.loc[end,'lon']
        lat2 = coor.loc[end,'lat']
        dis1,ind1 = misc_utils.near_euc(lfp_coor['lon'].values,lfp_coor['lat'].values,(lon1,lat1))
        dis2,ind2 = misc_utils.near_euc(lfp_coor['lon'].values,lfp_coor['lat'].values,(lon2,lat2))
        if (dis1 <= 0.01) & (dis2 <= 0.01):
            lfp_tree = lfp_tree.append(tree.loc[i,:])
    lfp_tree = lfp_tree[['link_no','start_pnt','end_pnt','frst_ds','frst_us','scnd_us','strahler','mon_pnt','shreve']]
    lfp_tree.index.name = 'index'

    # Writing clipped coord and tree files
    fnametre = folder + "/" + ncatchstr +"_tre.csv"
    fnamecoo = folder + "/" + ncatchstr +"_coo.csv"
    lfp_coor.to_csv(fnamecoo)
    lfp_tree.to_csv(fnametre,float_format='%i')

    # Creating XXX_rec.csv file
    fnamerec = folder + "/" + ncatchstr +"_rec.csv"
    connections(fnametre,fnamecoo,fnamerec)

    # Finding xmin, xmax, ymin, ymax based on river network points
    rec  = pd.read_csv(fnamerec)
    xmin = rec['lon'].min()
    xmax = rec['lon'].max()
    ymin = rec['lat'].min()
    ymax = rec['lat'].max()

    # Clipping rasters
    demarrcli,demgeocli = gdal_utils.clip_raster(demtif,xmin,ymin,xmax,ymax)
    accarrcli,accgeocli = gdal_utils.clip_raster(acctif,xmin,ymin,xmax,ymax)
    wtharrcli,wthgeocli = gdal_utils.clip_raster(wthtif,xmin,ymin,xmax,ymax)
    dirarrcli,dirgeocli = gdal_utils.clip_raster(dirtif,xmin,ymin,xmax,ymax)
    netarrcli,netgeocli = gdal_utils.clip_raster(nettif,xmin,ymin,xmax,ymax)
    catarrcli,catgeocli = gdal_utils.clip_raster(cattif,xmin,ymin,xmax,ymax)

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
    gdal_utils.writeRaster(demarrcli,fnamedem,demgeocli,"Float32",nodata)
    gdal_utils.writeRaster(accarrcli,fnameacc,accgeocli,"Float32",nodata)
    gdal_utils.writeRaster(netarrcli,fnamenet,netgeocli,"Float32",nodata)
    gdal_utils.writeRaster(wtharrcli,fnamewth,wthgeocli,"Float32",nodata)
    gdal_utils.writeRaster(dirarrcli,fnamedir,dirgeocli,"Float32",nodata)

def connections(treef,coorf,outfile):
    
    """
    Finds connections between links and sort them
    First finds the connections for all links, then sort them
    starting from the link with more downstream connections, later
    write separate files inclusing coordinates and links
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
    df_rec.to_csv(outfile)

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

    # Writing file
    df_rec.to_csv(outfile)

if __name__ == '__main__':
    split(sys.argv[1:])
