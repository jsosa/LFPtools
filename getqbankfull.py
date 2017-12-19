#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# date: 13/oct/2017
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

import sys,getopt,os,shutil,subprocess
import ConfigParser
import shapefile as sf
import pandas as pd
import jrc_utils
import matplotlib.pyplot as plt
import multiprocessing as mp
from gdal_utils import *
from pyproj import Proj, transform
from scipy import stats

import pdb
# pdb.set_trace()

def getqbankfull(argv):

    opts, args = getopt.getopt(argv,"i:")
    for o, a in opts:
        if o == "-i": inifile  = a

    config = ConfigParser.SafeConfigParser()
    config.read(inifile)
    
    netf   = str(config.get('getqbankfull','netf'))
    output = str(config.get('getqbankfull','output'))
    proj   = str(config.get('getqbankfull','proj'))
    aep    = np.float64(config.get('getqbankfull','aep'))
    thresh = np.float64(config.get('getqbankfull','thresh'))
    nproc  = np.float64(config.get('getqbankfull','nproc')) # number of cpus to use

    # coordinates for bank elevations are based in river network mask
    net   = get_gdal_data(netf)
    geo   = get_gdal_geo(netf)
    iy,ix = np.where(net>0)
    x     = geo[8][ix]
    y     = geo[9][iy]

    w = sf.Writer(sf.POINT)
    w.field('x')
    w.field('y')
    w.field('qbfull')

    # Split x and y in nproc parts
    split_x = np.array_split(x,nproc)
    split_y = np.array_split(y,nproc)

    # DEBUG DEBUG DEBUG
    # myq = calc_qbankfull(w,x,y,aep)
    # DEBUG DEBUG DEBUG

    # Define a queue
    queue = mp.Queue()

    # Setup a list of processes that we want to run
    processes = []
    processes = [mp.Process(target=calc_qbankfull_mp, args=(i,queue,w,split_x[i],split_y[i],aep)) for i in range(len(split_x))]

    # Run processes
    for p in processes:
        p.start()

    # Get process results from the queue
    results = [queue.get() for p in processes]
    
    # Retrieve results in a particular order
    results.sort()
    results = [r[1] for r in results]

    # Stack results horizontally
    myq = np.hstack(results)

    # Write qbankfull discharges
    for i in range(len(myq)):
        w.point(x[i],y[i])
        w.record(x[i],y[i],myq[i])

    w.save("%s.shp" % output)

    # write .prj file
    prj = open("%s.prj" % output, "w")
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj)
    prj.write(srs.ExportToWkt())
    prj.close()

    fmt      = "GTiff"
    nodata   = -9999
    bnkname1 = output+".shp"
    bnkname2 = output+".tif"
    subprocess.call(["gdal_rasterize","-a_nodata",str(nodata),"-of",fmt,"-tr",str(geo[6]),str(geo[7]),"-a","qbfull","-a_srs",proj,"-te",str(geo[0]),str(geo[1]),str(geo[2]),str(geo[3]),bnkname1,bnkname2])

def calc_qbankfull_mp(pos,queue,w,x,y,aep):

    myq = np.ones([len(x)])*-9999

    # iterate over every width x-y pair in the shapefile
    for i in range(len(x)):

        print "getqbankfull.py - " + str(len(x)-i)
        near_x,near_y = jrc_utils.find_nearest(x[i],y[i],1,10)

        if near_x != None:
            tserie = jrc_utils.get_data(near_x,near_y)
            myq[i] = qbfull_gev_mle(tserie,aep)
        else:
            print "NONE at " + str(x[i])+", "+str(y[i]) + "VERIFY!"
            myq[i] = np.nan

    queue.put((pos,myq))

def calc_qbankfull(w,x,y,aep):

    myq = np.ones([len(x)])*-9999

    # iterate over every width x-y pair in the shapefile
    for i in range(2700,2710):
    # for i in range(len(x)):

        print "getqbankfull.py - " + str(len(x)-i)
        near_x,near_y = jrc_utils.find_nearest(x[i],y[i],1,10)

        if near_x != None:
            tserie = jrc_utils.get_data(near_x,near_y)
            myq[i] = qbfull_gev_mle(tserie,aep)
            # print myq[i]
            # tserie.to_csv(str(i)+'.csv')
            # tserie.plot(); plt.show()
            print myq[i]

        else:
            print "NONE at " + str(x[i])+", "+str(y[i]) + "VERIFY!"
            myq[i] = np.nan

    return myq

def qbfull_mean_monthly(dsjrc):

    """
    Should be improved, currently not working!

    """

    dsjrc['mon'] = dsjrc.index.month
    monmean = dsjrc.groupby('mon').mean()
    
    stdup = monmean.values.mean() + monmean.values.std()
    stddw = monmean.values.mean() - monmean.values.std()
    idx   = np.where((monmean>=stddw) & (monmean<=stdup)) 

    score = []
    qu = monmean.values.mean()
    for i in idx[0]:
        qm = monmean.values[i][0]
        try:
            qml = monmean.values[i-1][0]
        except:
            qml = monmean.values[11][0]
        try:
            qmr = monmean.values[i+1][0]
        except:
            qmr = monmean.values[0][0]
        score.append((qu-0.5*qm-0.25*(qml+qmr),i))
    
    score.sort() # get the minimun score!
    imonth = score[0][1]
    month  = score[0][1] + 1 # index refer from zero to 11, then adding 1 indicates the month
    myq    = monmean.values[imonth][0]

    return myq

    # # DEBUG
    # monmean.plot()
    # plt.scatter(x=month,y=myq)
    # plt.show(block=False)
    # plt.pause(0.5)
    # plt.close()

def qbfull_gev_mle(tserie,aep):

    array = np.copy(tserie.discharge.values)
    array.sort() # array should be sorted

    # find parameters for GEV distribution by using MLE
    fit = stats.genextreme.fit(array)
    ppf = stats.genextreme.ppf(1-1/aep,*fit)

    # TRICK TO FIX NEGATIVE QBANKFULL DISCHARGE VALUES
    # just run again the subroutine but increasing 1% values in the time series
    if ppf<0:
        fit = stats.genextreme.fit(array*1.001)
        ppf = stats.genextreme.ppf(1-1/aep,*fit)

    myq = ppf

    return myq

def qbfull_gev_lmo(tseries,aep):

    """
    Should be improved, currently not working!
    
    """

    # # find parameters for GEV distribution by using L-Moments
    # xmom = lmom.samlmu(tserie,tserie.size,3)
    # para = lmom.pelgev(xmom)
    # print para

    return myq

if __name__ == '__main__':
    getqbankfull(sys.argv[1:])