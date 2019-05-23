#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: j.sosa@bristol.ac.uk / sosa.jeison@gmail.com

import os
import lfptools as lfp
from glob import glob
from shutil import copyfile
from subprocess import call

dem = './062/062_dem_lidar.tif'

# Calling lfp-getwidths
lfp.getwidths(thresh=0.02,
              output='./062/lfptools/062_wdt',
              recf='./062/062_rec.csv',
              netf='./062/062_net.tif',
              proj='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
              fwidth='./062/062_wth_grwl.tif')

# Calling bankelevs
lfp.getbankelevs(outlier='yes',
                 method='near',
                 hrnodata=-9999,
                 thresh=0.00416,
                 output='./062/lfptools/062_bnk',
                 netf='./062/062_net.tif',
                 hrdemf=dem,
                 proj='+proj=longlat + ellps=WGS84 + datum=WGS84 + no_defs')

# Calling fixelevs
lfp.fixelevs(method='yamazaki',
             source='./062/lfptools/062_bnk.shp',
             output='./062/lfptools/062_bnkfix',
             netf='./062/062_net.tif',
             recf='./062/062_rec.csv',
             proj='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

# Calling rasterreample
lfp.rasterresample(nproc=4,
                   outlier='yes',
                   method='mean',
                   hrnodata=-9999,
                   thresh=0.00416,
                   demf=dem,
                   netf='./062/062_net.tif',
                   output='./062/lfptools/062_dem30.tif')

# Calling getdepths
lfp.getdepths(proj='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
              netf='./062/062_net.tif',
              method='depth_geometry',
              output='./062/lfptools/062_dpt',
              wdtf='./062/lfptools/062_wdt.shp',
              r=0.12,
              p=0.78)

# Calling bedelevs
lfp.getbedelevs(bnkf='./062/lfptools/062_bnkfix.shp',
                dptf='./062/lfptools/062_dpt.shp',
                netf='./062/062_net.tif',
                output='./062/lfptools/062_bed',
                proj='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

# Creating a folder to save LISFLOOD-FP files
try:
    os.makedirs('./062/lisfloodfp/')
except FileExistsError:
    pass

# Calling buildmodel
lfp.buildmodel(parlfp='./062/lisfloodfp/062.par',
               bcilfp='./062/lisfloodfp/062.bci',
               bdylfp='./062/lisfloodfp/062.bdy',
               evaplfp='./062/lisfloodfp/062.evap',
               gaugelfp='./062/lisfloodfp/062.gauge',
               stagelfp='./062/lisfloodfp/062.stage',
               dembnktif='./062/lisfloodfp/062_dembnk.tif',
               dembnktif_1D='./062/lisfloodfp/062_dembnk_1D.tif',
               bedtif='./062/lfptools/062_bed.tif',
               wdttif='./062/lfptools/062_wdt.tif',
               runcsv='./062/062_dis.csv',
               demtif='./062/lfptools/062_dem30.tif',
               fixbnktif='./062/lfptools/062_bnkfix.tif',
               dirtif='./062/062_dir.tif',
               reccsv='./062/062_rec.csv',
               date1='1998-03-01',
               date2='1998-05-01')
[copyfile(i, './062/lisfloodfp/'+os.path.basename(i)) for i in glob('./062/lfptools/*.asc')]
