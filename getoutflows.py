#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# date: 22/apr/2017
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

import ConfigParser
import os,sys,getopt,shutil
import numpy as np
from osgeo import osr
import pyshp.shapefile as sf

from IPython.core.debugger import Pdb; pdb = Pdb()
# pdb.set_trace()

def getoutflows(argv):

    """
    This function uses the output from streamnet function from
    TauDEM, specifically the "coord" and "tree" files to adjust
    DEM values from rivers and tributaries for flood using the 
    algorithm bank4flood (1d)

    """

    opts, args = getopt.getopt(argv,"i:")
    for o, a in opts:
        if o == "-i": inifile  = a

    config = ConfigParser.SafeConfigParser()
    config.read(inifile)

    treef  = str(config.get('getoutflows','treef'))
    coordf = str(config.get('getoutflows','coordf'))
    shpd   = str(config.get('getoutflows','shpd'))
    proj   = str(config.get('getoutflows','proj'))

    coord  = np.genfromtxt(coordf, delimiter="\t")
    coord  = np.delete(coord,0,1) # remove first column, is an empty column
    tree   = np.genfromtxt(treef, delimiter="\t")
    tree   = np.delete(tree,0,1) # remove first column, is an empty column

    try: 
        os.makedirs(shpd)
    except OSError:
        if not os.path.isdir(shpd):
            raise
    fname = shpd+"/"+os.path.split(shpd)[1]
    
    w = sf.Writer(sf.POINT)
    w.field('outflow')

    count=0
    for i in range(tree.shape[0]):

        if tree[i,3] == -1: # is a downstream link
            end = int(tree[i,2])

            if coord[end,2] == 0: # if ending point, distance to outlet is zero

                count = count+1
                lon   = coord[end,0]
                lat   = coord[end,1]
                
                # write point in .shp file
                w.point(lon,lat)
                w.record(count)

    w.save("%s.shp" % fname)

    # write .prj file
    prj = open("%s.prj" % fname, "w")
    srs = osr.SpatialReference()
    srs.ImportFromProj4(proj)
    prj.write(srs.ExportToWkt())
    prj.close()

if __name__ == '__main__':

    getoutflows(sys.argv[1:])