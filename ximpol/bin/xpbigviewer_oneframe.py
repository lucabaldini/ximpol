#!/usr/bin/env python
from astropy.io import fits
import numpy
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.logging_ import logger, startmsg, abort
from ximpol.evt.binning import xBinnedCountSpectrum
from ximpol.evt.binning import xBinnedMap
from ximpol.evt.binning import xBinnedModulationCube
from ximpol.evt.subselect import xEventSelect
from ximpol.evt.binning import xEventBinningMCUBE
import os,sys
from matplotlib import rc
rc('text', usetex=True)

file_map='casa_cmap.fits'
evt_file_path='casa.fits'
file_selected_path='casa_sel.fits'
file_selected_cube_path='casa_sel_mcube.fits'
outfile=None
for i,a in enumerate(sys.argv):
    if '-rad'   in a: rad     = float(sys.argv[i+1])
    elif '-dec' in a: dec     = float(sys.argv[i+1])
    elif '-ra'  in a: ra      = float(sys.argv[i+1])
    elif '-o'   in a: outfile = sys.argv[i+1]
    pass

fig=xBinnedMap(file_map).plot(show=False,subplot=(1,2,1))#,figure=fig)
fig.show_circles(ra,dec,rad,lw=1)
rad*=60
evtSelect=xEventSelect(evt_file_path, ra=ra, dec=dec, rad=rad, outfile=file_selected_path,
                       emax=None, emin=None, mc=False, mcsrcid=[], phasemax=None, phasemin=None, phimax=None, phimin=None, tmax=None, tmin=None)
evtSelect.select()

evtBin=xEventBinningMCUBE(file_selected_path, ebins=1, outfile=file_selected_cube_path, evfile=file_selected_path,
                          emax=10.0,nypix=256, ebinfile=None, phasebins=50,ebinalg='LIN',
                          xref=None,phibins=75,nxpix=256,tbins=100,proj='TAN',tstart=None,tbinalg='LIN',algorithm='MCUBE',mc=False,binsz=2.5,yref=None,tbinfile=None,emin=1.0,tstop=None)
evtBin.bin_()

binModulation = xBinnedModulationCube(file_selected_cube_path)
binModulation.plot(show=False,xsubplot=1)
for fit in binModulation.fit_results:
    print fit
    angle      = fit.phase
    angle_err  = fit.phase_error

    visibility = fit.visibility

    scale=10.
    dx1=visibility/scale*numpy.cos(angle+10*angle_err)
    dy1=visibility/scale*numpy.sin(angle+10*angle_err)

    dx2=visibility/scale*numpy.cos(angle-10*angle_err)
    dy2=visibility/scale*numpy.sin(angle-10*angle_err)

    print fit.phase,angle,dx1,dy1,dx2,dy2

    fig.show_arrows(ra,dec,dx1,dy1,color='g',alpha=1,width=1)
    fig.show_arrows(ra,dec,-dx1,-dy1,color='g',alpha=1,width=1)

    fig.show_arrows(ra,dec,dx2,dy2,color='g',alpha=1,width=1)
    fig.show_arrows(ra,dec,-dx2,-dy2,color='g',alpha=1,width=1)
    pass

fig.save(outfile)
if outfile is None: plt.show()
