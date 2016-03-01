#!/usr/bin/env python
#
# Copyright (C) 2015--2016, the ximpol team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU GengReral Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


"""Module producing a polarization map for an extended source.
"""
import numpy
from astropy.io import fits
from astropy import wcs

from ximpol.utils.logging_ import logger, abort
#from ximpol.evt.event import xEventFile
#from ximpol.core.fitsio import xPrimaryHDU, xBinTableHDUBase
#from ximpol.irf.mrf import xAzimuthalResponseGenerator
from ximpol.utils.matplotlib_ import pyplot as plt

class xPolMapBase:
    ''' Base class for a polarization map
    '''
    def __init__(self, **kwargs):
        """Constructor.
        """
        self.kwargs = kwargs
        self.process_kwargs()
        
    def get(self, key, default=None):
        """Convenience method to address the keyword aguments.
        """
        return self.kwargs.get(key, default)
    
        
    def process_kwargs(self):
        pass

    def bin_(self):
        pass

    pass

class xPolMapRadial(xPolMapBase):
    ''' Radial polarization
    '''    
    def create(self):
        """Overloaded method.
        """
        xref = self.get('xref')
        yref = self.get('yref')
        nxpix = self.get('nxpix')
        pixsize = self.get('binsz')/3600.
        proj = self.get('proj')
        sidex = nxpix*pixsize
        sidey = nxpix*pixsize
        
        logger.info('Output image dimensions are %.1f x %.1f arcmin.' %\
                    (sidex*60, sidey*60))
                    
        logger.info('Center of the image is in R.A.=%.3f Dec.=%.3f' % (xref,yref))
                    
        binsx = numpy.linspace(0, nxpix, nxpix + 1)
        binsy = numpy.linspace(0, nxpix, nxpix + 1)
        # Build the WCS object
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = [nxpix/2., nxpix/2.]
        w.wcs.cdelt = [-pixsize, pixsize]
        w.wcs.crval = [xref, yref]

        #w.wcs.crpix = [0,0]
        #w.wcs.cdelt = [pixsize, pixsize]
        #w.wcs.crval = [xref - 0.5*sidex, yref - 0.5*sidey]
        w.wcs.ctype = ['RA---%s' % proj, 'DEC--%s' % proj]
        w.wcs.equinox = 2000.0
        #w.wcs.radesys = 'ICRS'
        header = w.to_header()
        
        pol_x=numpy.zeros((nxpix,nxpix))
        pol_y=numpy.zeros((nxpix,nxpix))
        
        for i in range(nxpix):
            for j in range(nxpix):
                #posx=w.wcs.crval[0]+i*pixsize
                #posy=w.wcs.crval[1]+j*pixsize
                # radial:
                #posx = nxpix/2.-i
                #posy = -nxpix/2.+j
                # circular:
                posx = -nxpix/2.+j
                posy = -nxpix/2.+i
                                
                ang = numpy.degrees(numpy.arctan2(posx,posy))
                deg = numpy.sqrt(posx*posx+posy*posy)
                                
                if (deg>0.5*nxpix) or i==0 or j==0 or i==(nxpix-1) or j==(nxpix-1):
                    posx=0
                    posy=0
                    pass
                pol_x[i,j]=posx
                pol_y[i,j]=posy           
                pass
            pass
        pol_x=numpy.array(pol_x)
        pol_y=numpy.array(pol_y)
        pol_deg=numpy.sqrt(pol_x*pol_x+pol_y*pol_y)
        pol_x/=pol_deg.max()
        pol_y/=pol_deg.max()
        outfile=self.get('outfile')
        outfile_x=outfile.replace('.fits','_x.fits')
        outfile_y=outfile.replace('.fits','_y.fits')
            
        hdu = fits.PrimaryHDU(pol_x, header=header)        
        logger.info('Writing binned Polarization X map in  %s...' % outfile_x)
        hdu.writeto(outfile_x, clobber=True)
        
        hdu = fits.PrimaryHDU(pol_y, header=header)        
        logger.info('Writing binned Polarization Y map in  %s...' % outfile_y)
        hdu.writeto(outfile_y, clobber=True)
        
        logger.info('Done.')
        pass
    pass

def readMap(filename):
    # Load the FITS hdulist using astropy.io.fits
    w = wcs.WCS(filename)
    
    # Print out the "name" of the WCS, as defined in the FITS header
    print(w.wcs.name)

    hdu = fits.open(filename)
    data=hdu[0].data
    nbinsx,nbinsy=data.shape
    pixcrd=[]
    data_arr=[]
    for x in range(nbinsx):
        for y in range(nbinsy):
            pixcrd.append([x,y])
            data_arr.append(data[x,y])
            pass
        pass
    world = w.wcs_pix2world(pixcrd, 0)
    return pixcrd, world, data_arr
    


if __name__ == '__main__':
    import os,sys
    from ximpol.srcmodel.img import xFITSImage
    from ximpol.utils.matplotlib_ import pyplot as plt
    from ximpol.utils.matplotlib_ import context_no_grids
    from ximpol.utils.logging_ import logger, startmsg
    outfile='out.fits'
    ra=None
    dec=None
    rad=None
    regionfile=None
    npix=20
    for i,a in enumerate(sys.argv):
        if '-o' in a: outfile    = sys.argv[i+1]
        if '-r' in a: regionfile = sys.argv[i+1]
        if '-ra' in a: ra=float(sys.argv[i+1])
        if '-dec' in a: dec=float(sys.argv[i+1])
        if '-rad' in a: rad=float(sys.argv[i+1])
        if '-npx' in a: npix= int(sys.argv[i+1])
        pass
    # center of the map, in RA, Dec:
    if regionfile is not None:
        for l in file(regionfile,'r').readlines():
            if 'circle(' in l:
                ra,dec,rad=l.replace('circle(','').replace('")','').split(',')
                break
            pass
        pass
    
    xref=float(ra)
    yref=float(dec)
    radius=float(rad)
    # Number of pixels:

    binsz  = 2.0*radius/npix
    
    myPolarizationMap = xPolMapRadial(xref=xref, yref=yref,nxpix=npix,nypix=npix,binsz=binsz,proj='TAN',outfile=outfile)
    myPolarizationMap.create()
    
    outfile_x=outfile.replace('.fits','_x.fits')
    outfile_y=outfile.replace('.fits','_y.fits')

    import aplpy
    import numpy
    import os
    from ximpol import XIMPOL_SRCMODEL

    #gc = aplpy.FITSFigure(outfile_deg,figsize=(10,9))
    gc = aplpy.FITSFigure(os.path.join(XIMPOL_SRCMODEL,'fits','casa_1p5_3p0_keV.fits'),figsize=(10,10))
    gc.show_colorscale(stretch='sqrt')#'graphics/2MASS_arcsinh_color.png')
    gc.set_tick_labels_font(size='small')
    #gc.show_contour(outfile_deg,colors='white')
    gc.show_grid()
    gc.add_scalebar(1./60.,'1\'')
    gc.scalebar.set_color('white')

    #gc.scalebar.show(0.2)  # length in degrees

    if regionfile is not None: gc.show_regions(regionfile)

    pixcrd, world, pol_x = readMap(outfile_x)
    pixcrd, world, pol_y = readMap(outfile_y)
    wx=[]
    wy=[]
    vx=[]
    vy=[]
    for i in range(world.shape[0]):
        x=world[i][0]
        y=world[i][1]
        px=0.01*pol_x[i]
        py=0.01*pol_y[i]
        if px*px+py*py>0:
            wx.append(x)
            wy.append(y)
            vx.append(px)
            vy.append(py)
            pass
        pass
    
    gc.show_arrows(wx,wy,vx,vy,color='w',alpha='0.8')
    gc.recenter(xref,yref,radius=1.2*radius/3600.0)

    plt.show()
