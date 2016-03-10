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
import aplpy
import os

from ximpol.utils.logging_ import logger, abort
#from ximpol.evt.event import xEventFile
#from ximpol.core.fitsio import xPrimaryHDU, xBinTableHDUBase
#from ximpol.irf.mrf import xAzimuthalResponseGenerator
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol import XIMPOL_CONFIG

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
    def add_circle(self,ra,dec,radius,pmax,type='radial'):
        
        for i in range(self.nxpix):
            for j in range(self.nxpix):
                # radial:
                #posx = nxpix/2.-i
                #posy = -nxpix/2.+j
                # circular:
                #p_x   = i #-nxpix+j
                #p_y   = j #-nxpix+i
                world = self.w.wcs_pix2world([[i,j]], 0)
                w_ra  = world[0][0]          
                w_dec = world[0][1]
                dx    =  (w_ra-ra)*numpy.cos(numpy.deg2rad(dec)) # Effect of the projection
                dy    =  w_dec-dec
                p_x = -dy
                p_y = +dx
                dist = numpy.sqrt(dx*dx+dy*dy)
                #print w_ra, w_dec, ra, dec, dist, radius
                if (dist>radius):
                    p_x=0
                    p_y=0
                    pass
                self.pol_x[i,j]=p_x
                self.pol_y[i,j]=p_y            
                pass
            pass
        pol_deg=numpy.sqrt(self.pol_x*self.pol_x+self.pol_y*self.pol_y)
        self.pol_x*=pmax/pol_deg.max()
        self.pol_y*=pmax/pol_deg.max()
        pass
    
    def save(self):
        outfile=self.get('outfile')
        outfile_x=outfile.replace('.fits','_x.fits')
        outfile_y=outfile.replace('.fits','_y.fits')            
        hdu = fits.PrimaryHDU(self.pol_x, header=self.header)        
        logger.info('Writing binned Polarization X map in  %s...' % outfile_x)
        hdu.writeto(outfile_x, clobber=True)
        hdu = fits.PrimaryHDU(self.pol_y, header=self.header)        
        logger.info('Writing binned Polarization Y map in  %s...' % outfile_y)
        hdu.writeto(outfile_y, clobber=True)
        logger.info('Done.')
        pass

    def create(self):
        """Overloaded method.
        """
        self.xref  = self.get('xref')
        self.yref  = self.get('yref')
        self.nxpix = self.get('nxpix')
        self.binsz = self.get('binsz') # Dimension of the pimap (degrees)
        self.proj = self.get('proj')

        
        sidex = self.nxpix*self.binsz
        logger.info('Output image dimensions are %.1f x %.1f arcmin.' %\
                    (sidex*60, sidex*60))
                    
        logger.info('Center of the image is in R.A.=%.3f Dec.=%.3f' % (self.xref,self.yref))                    
        # Build the WCS object
        self.w = wcs.WCS(naxis=2)
        self.w.wcs.crpix = [(self.nxpix+1)/2,(self.nxpix+1)/2]
        self.w.wcs.cdelt = [-self.binsz, self.binsz]
        self.w.wcs.crval = [self.xref, self.yref]
        self.w.wcs.ctype = ['RA---%s' % self.proj, 'DEC--%s' % self.proj]
        self.w.wcs.equinox = 2000.0
        #w.wcs.radesys = 'ICRS'
        self.header = self.w.to_header()
        self.pol_x=numpy.zeros((self.nxpix,self.nxpix))
        self.pol_y=numpy.zeros((self.nxpix,self.nxpix))
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
    import pyregion
    outfile=None
    ra=None
    dec=None
    rad=None
    regionfile=None
    npix=20
    pmax=1.0
    for i,a in enumerate(sys.argv):
        if '-o' in a: outfile    = sys.argv[i+1]
        if '-r' in a: regionfile = sys.argv[i+1]
        if '-ra' in a: ra=float(sys.argv[i+1])
        if '-dec' in a: dec=float(sys.argv[i+1])
        if '-rad' in a: rad=float(sys.argv[i+1])
        if '-npx' in a: npix= int(sys.argv[i+1])
        if '-pmax' in a: pmax= float(sys.argv[i+1])
        pass
    # First open the map (need for setting the geometry):
    le_img_file_path = os.path.join(XIMPOL_CONFIG, 'fits', 'casa_1p5_3p0_keV.fits')
    he_img_file_path = os.path.join(XIMPOL_CONFIG, 'fits', 'casa_4p0_6p0_keV.fits')
    img_file_path=le_img_file_path
    hdulist = fits.open(img_file_path)
    CRPIX1   = hdulist[0].header['CRPIX1']
    CRPIX2   = hdulist[0].header['CRPIX2']
    CDELT1   = hdulist[0].header['CDELT1']
    CDELT2   = hdulist[0].header['CDELT2']
    NAXIS1   = hdulist[0].header['NAXIS1']
    NAXIS2   = hdulist[0].header['NAXIS2']
    CRVAL1   = hdulist[0].header['CRVAL1']
    CRVAL2   = hdulist[0].header['CRVAL2']
    SIZE1    = CDELT1*NAXIS1
    SIZE2    = CDELT2*NAXIS2
    DEC_0    = CRVAL2+(NAXIS2/2-CRPIX2)*CDELT2
    # NOT SURE WHICH IS THE CORRECT PROJECTION...
    RA_0     = CRVAL1+(NAXIS1/2-CRPIX1)*CDELT1/numpy.cos(numpy.deg2rad(CRVAL2))
    #RA_0     = CRVAL1+(NAXIS1/2-CRPIX1)*CDELT1/numpy.cos(numpy.deg2rad(CRVAL2))

    gc = aplpy.FITSFigure(img_file_path,figsize=(10,10))
    gc.show_colorscale(stretch='sqrt')#'graphics/2MASS_arcsinh_color.png')
    gc.set_tick_labels_font(size='small')
    #gc.show_contour(outfile_deg,colors='white')
    gc.show_grid()
    gc.add_scalebar(1./60.,'1\'')
    gc.scalebar.set_color('white')
    gc.show_markers([CRVAL1],[CRVAL2],c='m')#,marker='x')
    gc.show_markers([RA_0],[DEC_0],c='g')#,marker='x')

    xref  = RA_0
    yref  = DEC_0
    binsz = max(SIZE1,SIZE2)/npix
    
    # center of the map, in RA, Dec:
    if regionfile is not None:
        if outfile is None:
            outfile = regionfile.replace('.reg','.fits')
        else:
            outfile='out.fits'
            pass
        regions = pyregion.open(regionfile)
        region = regions[0]
        ra, dec, rad = region.coord_list
        pass
    
    myPolarizationMap = xPolMapRadial(xref=xref, yref=yref,nxpix=npix,nypix=npix,binsz=binsz,proj='TAN',outfile=outfile,pmax=pmax)
    myPolarizationMap.create()
    myPolarizationMap.add_circle(ra,dec,rad,pmax)
    myPolarizationMap.save()
    #rad*=3600
    #logger.info(' Map center: %.3f,%.3f, radius=%.3f arcsec',ra,dec,rad)
    #xref=float(ra)
    #yref=float(dec)
    #radius=float(rad)
    # Number of pixels:
    #binsz  = 2.0*radius/npix
    
    
    outfile_x=outfile.replace('.fits','_x.fits')
    outfile_y=outfile.replace('.fits','_y.fits')


    #gc = aplpy.FITSFigure(outfile_deg,figsize=(10,9))
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
    gc.show_markers(wx,wy,c='r')
    gc.show_markers([xref],[yref],c='w')#,marker='x')

    #gc.recenter(xref,yref,radius=1.2*radius/3600.0)

    plt.show()
