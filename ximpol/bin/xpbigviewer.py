#!/usr/bin/env python
from astropy.io import fits
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.logging_ import logger, startmsg, abort
from ximpol.evt.binning import xBinnedCountSpectrum
from ximpol.evt.binning import xBinnedMap
from ximpol.evt.binning import xBinnedModulationCube
from ximpol.evt.select import xEventSelect
from ximpol.evt.binning import xEventBinningMCUBE
import os,sys
#from matplotlib import rc
#rc('text', usetex=True)

file_map='casa_cmap.fits'
evt_file_path='casa.fits'
import pyregion
region_name='casa_multiple.reg'
r = pyregion.open(region_name)
print r
number_of_regions=len(r)

#fig=plt.figure('my BIG canvas',figsize=(20,10))
#ax1=plt.subplot(1,2,1)
for i in range(number_of_regions):
    file_selected_path='casa_sel_%04d.fits' % i 
    file_selected_cube_path='casa_sel_mcube_%04d.fits' % i 
    if i<4: continue
    
    fig=xBinnedMap(file_map).plot(show=False,subplot=(1,2,1))#,figure=fig)
    ra,dec,rad=r[i].coord_list
    print '=====> ANALYZING FRAME %d: %f, %f, %f' %(i, ra,dec,rad)

    #ra=ras[i]
    #dec=decs[i]
    fig.show_circles(ra,dec,rad)
    rad*=60
    xEventSelect(evt_file_path, ra=ra, dec=dec, rad=rad, outfile=file_selected_path,
                 emax=None, emin=None, mc=False, mcsrcid=[], phasemax=None, phasemin=None, phimax=None, phimin=None, tmax=None, tmin=None).select()    

    xEventBinningMCUBE(file_selected_path, ebins=2, outfile=file_selected_cube_path,emax=10.0,nypix=256,evfile=file_selected_path, ebinfile=None, phasebins=50,ebinalg='LIN',
                       xref=None,phibins=75,nxpix=256,tbins=100,proj='TAN',tstart=None,tbinalg='LIN',algorithm='MCUBE',mc=False,binsz=2.5,yref=None,tbinfile=None,emin=1.0,tstop=None).bin_()
    xBinnedModulationCube(file_selected_cube_path).plot(xsubplot=1)
    fig.save('frame_%04d.png' %i)
    pass
plt.show()
