#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: j.sosa@bristol.ac.uk / sosa.jeison@gmail.com

import os
import sys
import getopt
import subprocess
import configparser
import numpy as np
import pandas as pd
import gdalutils


def buildmodel_shell(argv):
    """
    Main program to build a LISFLOOD-FP model
    """

    opts, args = getopt.getopt(argv, "i:")
    for o, a in opts:
        if o == "-i":
            inifile = a
    config = configparser.SafeConfigParser()
    config.read(inifile)

    runcsv = str(config.get('buildmodel', 'runcsv')) # input (discharge not from lfptools)

    demtif = str(config.get('buildmodel', 'demtif'))     # input 
    fixbnktif = str(config.get('buildmodel', 'fixbnktif')) # input 
    wdttif = str(config.get('buildmodel', 'wdttif')) # input 
    bedtif = str(config.get('buildmodel', 'bedtif')) # input 
    dirtif = str(config.get('buildmodel', 'dirtif')) # input 
    chantif = config.get('buildmodel', 'chantif',None) # input 
    reccsv = str(config.get('buildmodel', 'reccsv')) # input 
    date1 = str(config.get('buildmodel', 'date1')) # input 
    date2 = str(config.get('buildmodel', 'date2')) # input

    dembnktif = str(config.get('buildmodel', 'dembnktif')) # output 
    dembnktif_1D = str(config.get('buildmodel', 'dembnktif_1D')) # output 
    evaplfp = str(config.get('buildmodel', 'evaplfp')) # output
    gaugelfp = str(config.get('buildmodel', 'gaugelfp')) # output
    stagelfp = str(config.get('buildmodel', 'stagelfp')) # output
    parlfp = str(config.get('buildmodel', 'parlfp')) # output
    bcilfp = str(config.get('buildmodel', 'bcilfp')) # output 
    bdylfp = str(config.get('buildmodel', 'bdylfp')) # output
    d8dirn = config.getboolean('buildmodel','d8dirn',False)
    prescribeDirn = config.getboolean('buildmodel', 'prescribeDirn',False)

    buildmodel(parlfp, bcilfp, bdylfp, runcsv, evaplfp, gaugelfp, stagelfp,
               demtif, dembnktif, dembnktif_1D, fixbnktif, wdttif,
               bedtif, dirtif, reccsv, date1, date2, d8dir=d8dirn,prescribeDirn=prescribeDirn,chantif=chantif)


def buildmodel(parlfp, bcilfp, bdylfp, runcsv, evaplfp, gaugelfp, stagelfp,
               demtif, dembnktif, dembnktif_1D, fixbnktif, wdttif,
               bedtif, dirtif, reccsv, date1, date2, d8dirn=False,prescribeDirn=False,chantif=None):

    print("    running buildmodel.py...")

    t = (pd.to_datetime(date2, format='%Y-%m-%d') - pd.to_datetime(date1,
                                                                   format='%Y-%m-%d')).days + 1  # +1 to take into account the first date

    write_bci(bcilfp, runcsv)
    write_bdy(bdylfp, runcsv, t)
    write_evap(evaplfp, t)
    #write_gauge_stage_all_cells(reccsv, dirtif, wdttif, gaugelfp, stagelfp)
    burn_banks_dem(dembnktif, demtif, fixbnktif)
    burn_banks_dem_1D(dembnktif_1D, demtif, fixbnktif)
    write_ascii(dembnktif_1D, wdttif, bedtif, dembnktif,dirtif,chantif)
    if prescribeDirn:
        prescribeDirn = dirtif
    write_par(parlfp, bcilfp, bdylfp, evaplfp, gaugelfp,
              stagelfp, dembnktif, wdttif, bedtif, t,chantif,d8dirn,prescribeDirn)


def write_gauge_stage_all_cells(reccsv, dirtif, wdttif, gaugelfp, stagelfp):

    print("     writing gauge and stage files...")

    # Reading rec file
    rec = pd.read_csv(reccsv)

    # Create a width dataframe
    dat = gdalutils.get_data(wdttif)
    geo = gdalutils.get_geo(wdttif)
    wdt = gdalutils.array_to_pandas(dat, geo, 0, 'gt')
    wdt.columns = ['x', 'y', 'width']

    # Create directions dataframe
    dat = gdalutils.get_data(dirtif)
    geo = gdalutils.get_geo(dirtif)
    drc = gdalutils.array_to_pandas(dat, geo, 0, 'gt')
    drc.columns = ['x', 'y', 'direction']

    # Find widths and directions for every lon, lat in river network
    gdalutils.assign_val(df2=rec, df2_x='lon', df2_y='lat',
                         df1=wdt, df1_x='x', df1_y='y', label='width', copy=False)
    gdalutils.assign_val(df2=rec, df2_x='lon', df2_y='lat', df1=drc,
                         df1_x='x', df1_y='y', label='direction', copy=False)

    # Change numbers (1,2,3,4,5,6,7) to letters (N,S,E,W)
    rec['direction_let'] = rec['direction'].apply(getdirletter)

    # Writing .gauge file
    with open(gaugelfp, 'w') as f:
        f.write(str(rec.shape[0])+'\n')
    rec[['lon', 'lat', 'direction_let', 'width']].to_csv(
        gaugelfp, index=False, sep=' ', header=False, float_format='%.7f', mode='a')

    # Writing .stage file
    with open(stagelfp, 'w') as f:
        f.write(str(rec.shape[0])+'\n')
    rec[['lon', 'lat']].to_csv(
        stagelfp, index=False, sep=' ', header=False, float_format='%.7f', mode='a')


def write_evap(evaplfp, t):
    """
    writing Evaporation file
    Using 5mm/day evaporation value
    """

    print("     writing .evap file...")

    with open(evaplfp, "w") as f:
        f.write("# time series"+"\n")
        time = np.arange(t)  # daily values
        f.write(str(t)+"    "+"days"+"\n")
        for i in range(t):
            f.write("%12.3f    %d" % (10, time[i])+"\n")


def write_bdy(bdylfp, runcsv, t):
    """
    Subroutine used to write .BDY file
    Inflows are based on JRC hydrological model output
    """

    print("     writing .bdy file...")

    run = pd.read_csv(runcsv, index_col=0)

    # Select only date columns
    rund = run[[i for i in run.columns if (i[0] == '1') | (i[0] == '2')]].T

    # creating file
    with open(bdylfp, 'w') as f:
        f.write('# euflood bdy file'+'\n')

    # writing inflows
    for i in rund.columns:
        r = rund[i].to_frame()
        r['hours'] = range(0, t*24, 24)
        with open(bdylfp, 'a') as f:
            f.write('in'+str(i)+'\n')
            f.write(str(r['hours'].size)+' '+'hours'+'\n')
        r.to_csv(bdylfp, sep=' ', float_format='%.7f',
                 index=False, header=False, mode='a')


def write_bci(bcilfp, runcsv):
    """
    Writes bcif: XXX.bci file to be used in LISFLOOD-FP
    Uses runfcsv: XXX_run.csv
    """

    print("     writing .bci file...")

    run = pd.read_csv(runcsv, index_col=0)

    runi = run[['x', 'y']].T

    # creating file
    with open(bcilfp, 'w') as f:
        f.write('# euflood bci file'+'\n')

    # writing inflows
    with open(bcilfp, 'a') as f:
        for i in runi.columns:
            t = 'P'
            x = str(runi[i].loc['x'])
            y = str(runi[i].loc['y'])
            n = 'in' + str(i)
            f.write(t + ' ' + x + ' ' + y + ' ' + 'QVAR' + ' ' + n + '\n')
    
    # Writing other bundary conditions (eg. wall)

    with open(bcilfp,'a') as f:
        f.write('N -9999 9999 FREE' + '\n')
        f.write('S -9999 9999 FREE' + '\n')
        f.write('E -9999 9999 FREE' + '\n')
        f.write('W -9999 9999 FREE' + '\n')


def write_ascii(dembnktif_1D, wdttif, bedtif, dembnktif,dirtif,chantif):

    print("     writing ASCII files...")

    fmt2 = "AAIGRID"

    name2 = dembnktif_1D
    name3 = os.path.splitext(dembnktif_1D)[0]+'.asc'
    subprocess.call(["gdal_translate", "-of", fmt2,"-co","DECIMAL_PRECISION=2", name2, name3])

    name2 = wdttif
    name3 = os.path.splitext(wdttif)[0]+'.asc'
    subprocess.call(["gdal_translate", "-of", fmt2, name2, name3])

    name2 = bedtif
    name3 = os.path.splitext(bedtif)[0]+'.asc'
    subprocess.call(["gdal_translate", "-of", fmt2, name2, name3])

    name2 = dembnktif
    name3 = os.path.splitext(dembnktif)[0]+'.asc'
    subprocess.call(["gdal_translate", "-of", fmt2,"-co","DECIMAL_PRECISION=2", name2, name3])

    name2 = dirtif
    name3 = os.path.splitext(dirtif)[0]+'.asc'
    subprocess.call(["gdal_translate", "-of",fmt2,"-ot",'Int16', name2, name3])
	
    if chantif is not None:
        name2 = chantif
        name3 = os.path.splitext(chantif)[0]+'.asc'
        subprocess.call(["gdal_translate", "-of",fmt2,"-ot",'Int16', name2, name3])

def burn_banks_dem(dembnktif, demtif, fixbnktif):

    print("     burning banks in dem...")

    nodata = -9999
    fout = dembnktif
    base = gdalutils.get_data(demtif)
    basegeo = gdalutils.get_geo(demtif)
    new = gdalutils.get_data(fixbnktif)
    out = np.where(new > 0, new, base)
    gdalutils.write_raster(out, fout, basegeo, "Float32", nodata)


def burn_banks_dem_1D(dembnktif, demtif, fixbnktif):

    print("     burning banks in dem 1D...")

    nodata = -9999
    fout = dembnktif
    base = gdalutils.get_data(demtif)
    basegeo = gdalutils.get_geo(demtif)
    new = (np.ma.masked_values(gdalutils.get_data(
        fixbnktif), nodata)+10000).filled(nodata)
    out = np.where(new > 0, new, base)
    gdalutils.write_raster(out, fout, basegeo, "Float32", nodata)


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


def write_par(parlfp, bcilfp, bdylfp, evaplfp, gaugelfp, stagelfp, dembnktif, wdttif, bedtif, t, chantif,d8dirn,dirtif):

    print("     writing .par file...")

    with open(parlfp, "w") as file:

        file.write("latlong" + "\n")
        if d8dirn:
            file.write("SGCd8" + "\n")
        file.write("dirroot        " +
                   os.path.basename(parlfp).split('.')[0] + "\n")
        file.write("resroot        " +
                   os.path.basename(parlfp).split('.')[0] + "\n")
        file.write("sim_time       " + str((t-1)*86400) +
                   "\n")  # t-1, because first date
        file.write("initial_tstep  " + "10.0" + "\n")
        file.write("massint        " + "86400.0" + "\n")
        file.write("saveint        " + "86400.0" + "\n")
        file.write("fpfric         " + "0.06" + "\n")
        file.write("SGCn           " + "0.035" + "\n")
        if os.path.isfile(bcilfp):
            file.write("bcifile        " + os.path.basename(bcilfp) + "\n")
        if os.path.isfile(bdylfp):
            file.write("bdyfile        " + os.path.basename(bdylfp) + "\n")
        if os.path.isfile(evaplfp):
            file.write("evaporation    " + os.path.basename(evaplfp) + "\n")
		# PFU remove gauge and stage files
#        if os.path.isfile(gaugelfp):
#            file.write("gaugefile      " + os.path.basename(gaugelfp) + "\n")
#        if os.path.isfile(stagelfp):
#            file.write("stagefile      " + os.path.basename(stagelfp) + "\n")
        if os.path.isfile(dembnktif):
            file.write("DEMfile        " +
                       os.path.basename(dembnktif).split('.')[0] + '.asc' + "\n")
        if os.path.isfile(dembnktif):
            file.write("SGCbank        " +
                       os.path.basename(dembnktif).split('.')[0] + '.asc' + "\n")
        if os.path.isfile(wdttif):
            file.write("SGCwidth       " +
                       os.path.basename(wdttif).split('.')[0] + '.asc' + "\n")
        if os.path.isfile(bedtif):
            file.write("SGCbed         " +
                       os.path.basename(bedtif).split('.')[0] + '.asc' + "\n")
        if chantif is not None:
            if os.path.isfile(chantif):
                file.write("chanmask       " +
                       os.path.basename(chantif).split('.')[0] + '.asc' + "\n")
        if os.path.isfile(dirtif):
            file.write("SGCdirnfile    " +
                   os.path.basename(chantif).split('.')[0] + '.asc' + "\n")


if __name__ == '__main__':
    buildmodel_shell(sys.argv[1:])
