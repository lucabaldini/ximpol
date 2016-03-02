#!/urs/bin/env python
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

import numpy
from astropy.io import fits
from astropy import wcs
from ximpol.core.spline import xInterpolatedBivariateSplineLinear



def constant(C):
    """Simple wrapper returning a constant, independently of the input
    arguments.
    """
    def _function(E, t, ra, dec):
        return C
    return _function


class xPolMap:
    '''Class for the polarization map
    '''
    def __init__(self, x_filepath,y_filepath):
        self.x_wcs,self.x_map = self.readMap(x_filepath)
        self.y_wcs,self.y_map = self.readMap(y_filepath)

        pass

    def readMap(self,filename):
        # Load the FITS hdulist using astropy.io.fits
        w = wcs.WCS(filename)
        hdu = fits.open(filename)
        data=hdu[0].data
        nbinsx,nbinsy=data.shape
        x = range(nbinsx)
        y = range(nbinsy)
        spline=xInterpolatedBivariateSplineLinear(x,y,data)
        return w,spline

    def getPolarizationVector(self,ra,dec):
        """
        """
        x_x,x_y = self.x_wcs.wcs_world2pix(ra,dec,0)
        y_x,y_y = self.y_wcs.wcs_world2pix(ra,dec,0)
        px = self.x_map(x_x,x_y)
        py = self.y_map(y_x,y_y)
        return px,py

    def polarization_degree(self, ra, dec):
        """Return the polarization degree for a given direction in the sky.

        Warning
        -------
        Note that we're calling the getPolarizationVector() function here and
        in the polarization_angle() method, while we could in principle
        get away with just one function call. The issue is that downstream
        we need the polarization degree and angle separately, and if we want
        to optimize things here, we would have to implement a generic
        polarization() interface in the model components that is called
        consistently in rvs_event_list().
        """
        px, py = self.getPolarizationVector(ra, dec)
        return numpy.sqrt(px*px + py*py)

    def polarization_angle(self, ra, dec):
        """Return the polarization angle for a given direction in the sky.
        """
        px, py = self.getPolarizationVector(ra, dec)
        phi = numpy.arctan2(py, px)
        phi += (phi < 0.)*numpy.pi
        return phi


if __name__=='__main__':
    import os
    import aplpy
    from ximpol import XIMPOL_SRCMODEL
    from ximpol.utils.matplotlib_ import pyplot as plt
    import numpy.random

    polarization_x_map = os.path.join(XIMPOL_SRCMODEL,'fits','casa_pol_x.fits')
    polarization_y_map = os.path.join(XIMPOL_SRCMODEL,'fits','casa_pol_y.fits')
    casa_map           = os.path.join(XIMPOL_SRCMODEL,'fits','casa_1p5_3p0_keV.fits')
    polarization_map   = xPolMap(polarization_x_map,polarization_y_map)
    my_ra  = 350.863
    my_dec =  58.815
    npoints=1000
    delta=0.05
    #ras  = numpy.linspace(my_ra-delta,my_ra+delta,npoints)
    #decs = numpy.linspace(my_dec-delta,my_dec+delta,npoints)
    ras   = numpy.random.uniform(my_ra -2*delta,my_ra+2*delta,npoints)
    decs  = numpy.random.uniform(my_dec-delta,my_dec+delta,npoints)
    wx=[]
    wy=[]
    vx=[]
    vy=[]
    for ra,dec in zip(ras,decs):
        px,py = polarization_map.getPolarizationVector(ra,dec)
        if (px*px+py*py)>0:
            wx.append(ra)
            wy.append(dec)
            vx.append(0.01*px)
            vy.append(0.01*py)
            pass
        pass
    #aplpy.make_rgb_cube([polarization_x_map,polarization_y_map,casa_map],'thecube.fits')
    #aplpy.make_rgb_image('thecube.fits','test.tif',stretch_r='linear',stretch_g='linear',stretch_b='arcsinh')
    #aplpy.show_arrows(wx,wy,vx,vy,color='w',alpha=0.8)

    gc = aplpy.FITSFigure(polarization_x_map,figsize=(10,10))
    gc.show_colorscale(stretch='sqrt')
    gc.show_arrows(wx,wy,vx,vy,color='w',alpha=0.8)
    gc.show_markers(wx,wy,marker='o')
    #gc.tick_labels.set_xformat('dd')
    #gc.tick_labels.set_yformat('dd')

    gc = aplpy.FITSFigure(polarization_y_map,figsize=(10,10))
    gc.show_colorscale(stretch='sqrt')
    gc.show_arrows(wx,wy,vx,vy,color='w',alpha=0.8)
    gc.show_markers(wx,wy,marker='o')

    gc = aplpy.FITSFigure(casa_map,figsize=(10,10))
    gc.show_colorscale(stretch='sqrt')
    gc.show_arrows(wx,wy,vx,vy,color='w',alpha=0.8)
    gc.show_markers(wx,wy,marker='o')
    #gc.tick_labels.set_xformat('dd')
    #gc.tick_labels.set_yformat('dd')

    plt.show()
