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

class xPolMap(xPolMapBase):
    ''' Radial polarization
    '''
    def add_circle(self,ra,dec,radius,pmax,ptype='circular',angle=0.0):
        ''' This implement the case of a circular shape'''
        logger.info('add a circle at ra=%f, dec=%f, radius=%f. PMAX=%f'%(ra,dec,rad,pmax))
        for i in range(self.nxpix):
            for j in range(self.nxpix):
                world = self.w.wcs_pix2world([[i,j]], 0)
                w_ra  = world[0][0]          
                w_dec = world[0][1]
                dx    =  (w_ra-ra)*numpy.cos(numpy.deg2rad(dec)) # Effect of the projection
                dy    =  w_dec-dec
                if ptype is 'linear':
                    p_x = numpy.cos(numpy.deg2rad(angle))
                    p_y = numpy.sin(numpy.deg2rad(angle))
                elif ptype is 'circular':
                    p_x = -dy
                    p_y = +dx
                elif ptype is 'radial':
                    p_x = dx
                    p_y = dy
                    pass
                dist = numpy.sqrt(dx*dx+dy*dy)
                #print w_ra, w_dec, ra, dec, dist, radius
                if (dist>radius):
                    p_x=0
                    p_y=0
                    pass
                self.pol_x[i,j]=self.pol_x[i,j]+p_x
                self.pol_y[i,j]=self.pol_y[i,j]+p_y            
                pass
            pass
        pol_deg=numpy.sqrt(self.pol_x*self.pol_x+self.pol_y*self.pol_y)
        self.pol_x*=pmax/pol_deg.max()
        self.pol_y*=pmax/pol_deg.max()
        pass
    
    def add_ellipse(self,ra,dec,radius1,radius2,theta, pmax, ptype='circular',angle=0.0):
        ''' This implement the case of a circular shape'''
        logger.info('add a ellipse at ra=%f, dec=%f, radius1=%f radius2=%f theta=%f. PMAX=%f, ANGLE=%f, PTYPE=%s.'%(ra,dec,radius1,radius2,theta, pmax,angle,ptype))
        theta_rad=numpy.deg2rad(theta)
        for i in range(self.nxpix):
            for j in range(self.nxpix):
                world = self.w.wcs_pix2world([[i,j]], 0)
                w_ra  = world[0][0]
                w_dec = world[0][1]
                dx0   = (w_ra-ra)*numpy.cos(numpy.deg2rad(dec))
                dy0   = (w_dec-dec)
                #dx1   = dx0#dx0*numpy.cos(theta_rad)+dy0*numpy.sin(theta_rad)
                #dy1   = dy0#-dx0*numpy.sin(theta_rad)+dy0*numpy.cos(theta_rad)
                dx1   =  dx0*numpy.cos(theta_rad)-dy0*numpy.sin(theta_rad)
                dy1   =  dx0*numpy.sin(theta_rad)+dy0*numpy.cos(theta_rad)
                dx    =  dx1/radius1 # Effect of the projection
                dy    =  dy1/radius2
                dx2 = dx*numpy.cos(theta_rad)+dy*numpy.sin(theta_rad)
                dy2 = -dx*numpy.sin(theta_rad)+dy*numpy.cos(theta_rad)
                
                if ptype == 'linear':
                    p_x = numpy.cos(numpy.deg2rad(angle))
                    p_y = numpy.sin(numpy.deg2rad(angle))
                elif ptype == 'circular':
                    p_x = -dy2
                    p_y = +dx2
                elif ptype == 'radial':
                    p_x = dx2
                    p_y = dy2
                    pass
                dist = numpy.sqrt(dx*dx+dy*dy)
                if (dist>1):
                    p_x=0
                    p_y=0
                    pass
                else:
                    #print w_ra, w_dec, ra, dec, dist, dx1,dy1, dx2,dy2,p_x, p_y
                    pass
                self.pol_x[i,j]=self.pol_x[i,j]+p_x#*radius1
                self.pol_y[i,j]=self.pol_y[i,j]+p_y#*radius2           
                pass
            pass
        pol_deg=numpy.sqrt(self.pol_x*self.pol_x+self.pol_y*self.pol_y)
        self.pol_x*=pmax/pol_deg.max()
        self.pol_y*=pmax/pol_deg.max()
        pass
    
    def add_shape(self,myfilter,pmax,angle=0.0):
        ''' This is a generic shape (box, poliline' used to define a region of uniform polarization'''
        logger.info('add a spatial shape')
        for i in range(self.nxpix):
            for j in range(self.nxpix):
                #print i,j,myfilter.inside1(i, j)
                if myfilter.inside1(i, j): 
                    p_x = pmax*numpy.cos(numpy.deg2rad(angle))
                    p_y = pmax*numpy.sin(numpy.deg2rad(angle))
                    self.pol_x[i,j]=self.pol_x[i,j]+p_x
                    self.pol_y[i,j]=self.pol_y[i,j]+p_y  
                    pass
                pass
            pass
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
    import argparse
    desc = '''Toll for creating polarization maps from DS9 Region files'''
    parser = argparse.ArgumentParser(description=desc)
    
    parser.add_argument('--image',help='Input fits image file', type=str, required=True)
    parser.add_argument('--region',help='Input Region file containing DS9 regions', type=str, required=True)
    parser.add_argument('--output',help='Name for the output file (for the polarizationmap[s])',required=True, type=str)
    parser.add_argument('--npix',help='Number of points in the final polarization map',type=int, required=False, default=50)


    args = parser.parse_args()
        
    out_file_path=args.output
    img_file_path=args.image
    reg_file_path=args.region
    npix=args.npix
    # I need to find a good way to pass this value for different regions
    # pmax=0.5
    ##############################
    # First open the map (need for setting the geometry):
    logger.info('open image file: %s' % img_file_path)
    hdulist = fits.open(img_file_path)
    CRPIX1   = hdulist[0].header['CRPIX1']
    CRPIX2   = hdulist[0].header['CRPIX2']
    CDELT1   = hdulist[0].header['CDELT1']
    CDELT2   = hdulist[0].header['CDELT2']
    NAXIS1   = hdulist[0].header['NAXIS1']
    NAXIS2   = hdulist[0].header['NAXIS2']
    CRVAL1   = hdulist[0].header['CRVAL1']
    CRVAL2   = hdulist[0].header['CRVAL2']
    SIZE1    = numpy.fabs(CDELT1*NAXIS1)
    SIZE2    = numpy.fabs(CDELT2*NAXIS2)
    DEC_0    = CRVAL2+(NAXIS2/2-CRPIX2)*CDELT2
    #  PROJECTION...
    RA_0     = CRVAL1+(NAXIS1/2-CRPIX1)*CDELT1/numpy.cos(numpy.deg2rad(CRVAL2))

    gc = aplpy.FITSFigure(img_file_path,figsize=(10,10))
    gc.show_colorscale(stretch='sqrt')#'graphics/2MASS_arcsinh_color.png')
    gc.set_tick_labels_font(size='small')
    #gc.show_contour(outfile_deg,colors='white')
    gc.show_grid()
    gc.add_scalebar(1./60.,'1\'')
    gc.scalebar.set_color('white')
    #gc.show_markers([CRVAL1],[CRVAL2],c='m')#,marker='x')
    #gc.show_markers([RA_0],[DEC_0],c='g')#,marker='x')
    xref  = RA_0
    yref  = DEC_0
    binsz = max(SIZE1,SIZE2)/npix

    #if ra is None:  ra=RA_0
    #if dec is None: dec=DEC_0
    #if rad is None: rad=min(SIZE1,SIZE2)

    # center of the map, in RA, Dec:
    #if regionfile is not None:
    #    if outfile is None:
    #        outfile = regionfile.replace('.reg','_%d.fits' %(100*pmax))
    #    else:
    #        outfile='out_%d.fits' %(100*pmax)
    #        pass
    #pass

    # REGION FILE:
    logger.info('open image file: %s' % reg_file_path)
    # Open the region file
    myPolarizationMap = xPolMap(xref=xref, yref=yref,nxpix=npix,nypix=npix,binsz=binsz,proj='TAN')
    myPolarizationMap.create()
    regions_img   = pyregion.open(reg_file_path).as_imagecoord(myPolarizationMap.header)
    regions_word  = pyregion.open(reg_file_path)
    myfilters = regions_img.get_filter()
    Nregion = len(regions_img)
    out_files=[]
    for i in range(Nregion):
        r=regions_word[i]
        # I need to create one polarization map per region
        print '===> region found... NAME:', r.name,'COORD SYS:',r.coord_format,'COMMENT:',r.comment
        print '     COORD LIST:',r.coord_list
        if r.comment is not None and 'linear' in r.comment:
            ptype='linear'
        elif r.comment is not None and 'circular' in r.comment:
            ptype='circular'
        elif r.comment is not None and 'radial' in r.comment:
            ptype='radial'
            pass
        else:
            ptype=raw_input('===> Type of polarization pattern [linear|circular|radial]....: ').replace('\n','').strip()
            pass
        pangle=0
        if ptype=='linear':
            if r.comment is not None and 'pangle=' in r.comment:
                pangle=float(r.comment.split('pangle=')[-1].split(' ')[0])
                print 'Angle of polarization = %f degrees' % pangle
            else:
                pangle=float(raw_input('===> Angle in degrees....: '))    
                pass
            pass
        
        if r.comment is not None and 'pmax=' in r.comment:
            pmax=float(r.comment.split('pmax=')[-1].split(' ')[0])
            print 'Maximum Degree of polarization of polarization = %f' % pmax
        else:
            pmax=float(raw_input('===> maximum degree of polarization for region %s [0-100]....: ' %  r.name))/100.        
            pass

        outfile= out_file_path.replace('.fits','_pmax%03d_reg%03d.fits' % (pmax*100,i))
        out_files.append(outfile)
        myPolarizationMap = xPolMap(xref=xref, yref=yref,nxpix=npix,nypix=npix,binsz=binsz,proj='TAN',outfile=outfile)
        myPolarizationMap.create()
        #pmax*=Nregion # This is to account that in the simulation each 
        if r.name is 'circle':
            ra, dec, rad = r.coord_list
            myPolarizationMap.add_ellipse(ra,dec,rad,rad,0.0,pmax,ptype, pangle)
            #myPolarizationMap.add_circle(ra,dec,rad,pmax,ptype)
            pass
        elif r.name is 'ellipse':
            ra, dec, rad1, rad2, ang = r.coord_list
            myPolarizationMap.add_ellipse(ra,dec,rad1,rad2,ang,pmax,ptype, pangle)            
        else:
            myPolarizationMap.add_shape(myfilters[i],pmax=pmax,angle=pangle)
            pass
        myPolarizationMap.save()
        pass
    
    #Load THE OUTPUT:
    # DISPLAT THE REGION FILE
    gc.show_regions(reg_file_path)    
    for j,outfile in enumerate(out_files):
        outfile_x= outfile.replace('.fits','_x.fits') 
        outfile_y= outfile.replace('.fits','_y.fits') 
        #gc = aplpy.FITSFigure(outfile_deg,figsize=(10,9))
        #gc.scalebar.show(0.2)  # length in degrees
        # DIPLAY THE POLARIZATION MAP:
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
        colors=['w','y','g','r','b','k']
        _=gc.show_arrows(wx,wy,vx,vy,color=colors[j % len(colors)],alpha='1.0',linewidth=2)
        #gc.show_markers(wx,wy,c='r')
        #gc.show_markers([xref],[yref],c='w')#,marker='x')
        #gc.recenter(xref,yref,radius=1.2*radius/3600.0)
        pass
    plt.show()
