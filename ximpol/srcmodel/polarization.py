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
from ximpol.srcmodel.img import xFITSImage
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.matplotlib_ import context_no_grids

def constant(C):
    """Simple wrapper returning a constant, independently of the input
    arguments.
    """
    def _function(E, t, ra, dec):
        return C
    return _function


class xPolarizationMap:
    """Read-mode interface for the polarization maps.
    """

    def __init__(self, xmap_file_path, ymap_file_path):
        """Constructor.
        """
        self.xmap_file_path = xmap_file_path
        self.ymap_file_path = ymap_file_path
        self.x_img = None
        self.y_img = None
        self.x_wcs, self.x_data, self.x_spline = self.__read_map(xmap_file_path)
        self.y_wcs, self.y_data, self.y_spline = self.__read_map(ymap_file_path)
        self.__reset_sample()

    def __read_map(self,filename):
        """Load the FITS hdulist using astropy.io.fits.
        """
        w = wcs.WCS(filename)
        hdu = fits.open(filename)
        data = hdu[0].data
        nbinsx, nbinsy = data.shape
        x = range(nbinsx)
        y = range(nbinsy)
        spline = xInterpolatedBivariateSplineLinear(x, y, data)
        return w, data, spline

    def polarization_vector(self, ra, dec):
        """Return the polarization vector for a given (array of) RA and Dec
        values.
        """
        x_x, x_y = self.x_wcs.wcs_world2pix(ra, dec, 0)
        y_x, y_y = self.y_wcs.wcs_world2pix(ra, dec, 0)
        px = self.x_spline(x_x, x_y)
        py = self.y_spline(y_x, y_y)
        return px, py

    def polarization_degree(self, ra, dec):
        """Return the polarization degree for a given direction in the sky.

        Warning
        -------
        Note that we're calling the polarization_vector() function here and
        in the polarization_angle() method, while we could in principle
        get away with just one function call. The issue is that downstream
        we need the polarization degree and angle separately, and if we want
        to optimize things here, we would have to implement a generic
        polarization() interface in the model components that is called
        consistently in rvs_event_list().
        """
        px, py = self.polarization_vector(ra, dec)
        return numpy.sqrt(px*px + py*py)

    def polarization_angle(self, ra, dec):
        """Return the polarization angle for a given direction in the sky.
        """
        px, py = self.polarization_vector(ra, dec)
        phi = numpy.arctan2(py, px)
        phi += (phi < 0.)*numpy.pi
        return phi

    def __reset_sample(self):
        """
        """
        self.__wx = []
        self.__wy = []
        self.__vx = []
        self.__vy = []

    def build_random_sample(self, ra0, dec0, num_points=1000, radius=5.):
        """Calculate the polarization components on a random set of
        points.

        Warning
        -------
        This could be probably vectorized.
        """
        self.__reset_sample()
        radius /= 60.
        delta_ra = radius/numpy.cos(numpy.radians(dec0))
        ra = numpy.random.uniform(ra0 - delta_ra, ra0 + delta_ra, num_points)
        dec = numpy.random.uniform(dec0 - radius, dec0 + radius, num_points)
        for _ra, _dec in zip(ra, dec):
            px, py = self.polarization_vector(_ra, _dec)
            if (px*px + py*py) > 0:
                self.__wx.append(_ra)
                self.__wy.append(_dec)
                self.__vx.append(0.01*px)
                self.__vy.append(0.01*py)

    def build_grid_sample(self, ra0, dec0, num_points=25, radius=5.):
        """Calculate the polarization components on a rectangular grid.

        Warning
        -------
        This could be probably vectorized.
        """
        self.__reset_sample()
        radius /= 60.
        delta_ra = radius/numpy.cos(numpy.radians(dec0))
        ra = numpy.linspace(ra0 - delta_ra, ra0 + delta_ra, num_points)
        dec = numpy.linspace(dec0 - radius, dec0 + radius, num_points)
        for _ra in ra:
            for _dec in dec:
                px, py = self.polarization_vector(_ra, _dec)
                if (px*px + py*py) > 0:
                    self.__wx.append(_ra)
                    self.__wy.append(_dec)
                    self.__vx.append(0.01*px)
                    self.__vy.append(0.01*py)

    def overlay_arrows(self, fig, markers=True):
        """Overlay the polarization map arrows over an existing aplpy figure.
        """
        fig.show_arrows(self.__wx, self.__wy, self.__vx, self.__vy,
                        color='w', alpha=0.8)
        if markers:
            fig.show_markers(self.__wx, self.__wy, marker='o')

    def plot_xmap(self, overlay=True, show=False):
        """Plot the x polarization map.
        """
        if self.x_img is None:
            self.x_img = xFITSImage(self.xmap_file_path, build_cdf=False)
        fig = self.x_img.plot(show=False, zlabel='Polarization degree (x)')
        if overlay:
            self.overlay_arrows(fig)
        if show:
            plt.show()
        return fig

    def plot_ymap(self, overlay=True, show=False):
        """Plot the y polarization map.
        """
        if self.y_img is None:
            self.y_img = xFITSImage(self.ymap_file_path, build_cdf=False)
        fig = self.y_img.plot(show=False, zlabel='Polarization degree (y)')
        if overlay:
            self.overlay_arrows(fig)
        if show:
            plt.show()
        return fig

    def plot_polarization_degree(self, show=True):
        """
        """
        import aplpy
        if self.x_img is None:
            self.x_img = xFITSImage(self.xmap_file_path, build_cdf=False)
        if self.y_img is None:
            self.y_img = xFITSImage(self.ymap_file_path, build_cdf=False)
        _data = numpy.sqrt(self.x_data**2 + self.y_data**2)
        hdu_list = [self.x_img.hdu_list[0].copy()]
        hdu_list[0].data = _data
        with context_no_grids():
            fig = aplpy.FITSFigure(hdu_list[0], figure=plt.figure())
            fig.add_grid()
            fig.show_colorscale(cmap = 'afmhot', vmin=None, vmax=None)
            fig.add_colorbar()
            fig.colorbar.set_axis_label_text('Polarization degree')
        if show:
            plt.show()
        return fig
        
    def plot_polarization_angle(self, degrees=True, show=True):
        """
        """
        import aplpy
        if self.x_img is None:
            self.x_img = xFITSImage(self.xmap_file_path, build_cdf=False)
        if self.y_img is None:
            self.y_img = xFITSImage(self.ymap_file_path, build_cdf=False)
        _data = numpy.arctan2(self.x_data, self.y_data)
        if degrees:
            _data = numpy.degrees(_data)
        hdu_list = [self.x_img.hdu_list[0].copy()]
        hdu_list[0].data = _data
        with context_no_grids():
            fig = aplpy.FITSFigure(hdu_list[0], figure=plt.figure())
            fig.add_grid()
            fig.show_colorscale(cmap = 'afmhot', vmin=None, vmax=None)
            fig.add_colorbar()
            fig.colorbar.set_axis_label_text('Polarization angle')
        if show:
            plt.show()
        return fig
    


if __name__=='__main__':
    import os
    from ximpol import XIMPOL_CONFIG
    file_path_x = os.path.join(XIMPOL_CONFIG, 'fits', 'casa_pol_x.fits')
    file_path_y = os.path.join(XIMPOL_CONFIG, 'fits', 'casa_pol_y.fits')
    img_file_path = os.path.join(XIMPOL_CONFIG, 'fits', 'casa_1p5_3p0_keV.fits')
    polarization_map = xPolarizationMap(file_path_x, file_path_y)
    polarization_map.build_grid_sample(350.863, 58.815)
    polarization_map.plot_xmap()
    polarization_map.plot_ymap()
    img = xFITSImage(img_file_path)
    fig = img.plot(show=False)
    polarization_map.overlay_arrows(fig)
    plt.show()
