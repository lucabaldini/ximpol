#!/usr/bin/env python
import os,sys
import pyregion
region_name=sys.argv[1]
r = pyregion.open(region_name)
number_of_regions=len(r)
print 'found %d regions...' % number_of_regions
for i in range(number_of_regions):
    ra,dec,rad=r[i].coord_list
    outfilename='frame_%04d.png' %i
    cmd='./xpbigviewer_oneframe.py -ra %s -dec %s -rad %s -o %s' %(ra,dec,rad,outfilename)
    print cmd
    os.system(cmd)
    pass
