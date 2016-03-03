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


from astropy.io import fits
from astropy.wcs import wcs
import numpy

from ximpol.utils.logging_ import logger
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.matplotlib_ import context_no_grids


class xFITSImage:

    """Class describing a FITS image.

    Arguments
    ---------
    file_path : string
        The path to the FITS file containing the image.

    build_cdf : bool
        If True, build the cdf (i.e., equip the instance to generate random
        numbers).

    Warning
    -------
    There are several things I don't quite understand here, first of all
    why we seem to need to transpose the data. (Also, we migh have a
    residual offset by 1 pixel that we should try and sort out.)
    """

    def __init__(self, file_path, build_cdf=True):
        """Constructor.
        """
        logger.info('Reading FITS image from %s...' % file_path)
        self.hdu_list = fits.open(file_path)
        self.hdu_list.info()
        self.wcs = wcs.WCS(self.hdu_list['PRIMARY'].header)
        self.data = self.hdu_list['PRIMARY'].data.transpose()
        if build_cdf:
            self.build_cdf()

    def build_cdf(self):
        """Build the cumulative distribution function.

        (This is used to extract random positions from the image when
        simulating extended sources.)
        """
        self.cdf = numpy.cumsum(self.data.ravel())
        self.cdf /= self.cdf[-1]

    def rvs_coordinates(self, size=1, randomize=True):
        """Generate random coordinates based on the image map.

        Arguments
        ---------
        size : int
            The number of sky coordinates to be generated.

        randomize : bool
            If true, the positions are randomized uniformely within each pixel.

        Warning
        -------
        There must be a better way to do this. We should take a look at
        how aplpy.FITSImage is doing this.
        """
        u = numpy.random.rand(size)
        pixel = numpy.searchsorted(self.cdf, u)
        row, col = numpy.unravel_index(pixel, self.data.shape)
        pixel_crd = numpy.vstack((row, col)).transpose()
        world_crd = self.wcs.wcs_pix2world(pixel_crd, 1)
        ra, dec = world_crd[:, 0], world_crd[:, 1]
        if randomize:
            delta_ra = 0.5*self.hdu_list['PRIMARY'].header['CDELT1']
            delta_dec = 0.5*self.hdu_list['PRIMARY'].header['CDELT2']
            ra += numpy.random.uniform(-delta_ra, delta_ra, size)
            dec += numpy.random.uniform(-delta_dec, delta_dec, size)
        return ra, dec

    def __call__(self, row, column):
        """Return the value of the underlying map for a given pixel.
        """
        return self.data[i][j]

    def plot(self, show=True, zlabel='Counts/pixel',subplot=(1,1,1)):
        """Plot the image.

        This is using aplpy to render the image.
        """
        import aplpy
        with context_no_grids():
            fig = aplpy.FITSFigure(self.hdu_list[0], figure=plt.figure(0,figsize=(10*subplot[1],10*subplot[0])), subplot = subplot)
        fig.add_grid()
        fig.show_colorscale(cmap = 'afmhot')
        fig.add_colorbar()
        fig.colorbar.set_axis_label_text(zlabel)
        if show:
            plt.show()
        return fig



def main():
    """
    """
    import os
    from ximpol import XIMPOL_SRCMODEL
    file_path = os.path.join(XIMPOL_SRCMODEL, 'fits', 'crab_0p3_10p0_keV.fits')
    img = xFITSImage(file_path)
    ra, dec = img.rvs_coordinates(1000000)
    print(ra)
    print(dec)
    img.plot()


if __name__ == '__main__':
    main()
