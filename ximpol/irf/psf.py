#!/usr/bin/env python
#
# Copyright (C) 2015, the ximpol team.
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

from ximpol.utils.logging_ import logger
from ximpol.irf.base import xColDefsBase, OGIP_HEADER_SPECS
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear


PSF_HEADER_SPECS = [
    ('EXTNAME' , 'PSF', 'name of this binary table extension')
]


class xColDefsPSF(xColDefsBase):

    """ximpol.irf.base.xColDefsBase subclass for the PSF extension
    of .psf FITS files.
    """

    COLUMN_SPECS = [
        ('W'    , 'E', '1/sr'),
        ('SIGMA', 'E', 'arcsec'),
        ('N'    , 'E', '1/sr'),
        ('R_C'  , 'E', 'arcsec'),
        ('ETA'  , 'E', None)
    ]


class xPointSpreadFunction(xInterpolatedUnivariateSplineLinear):

    """Class describing a (simplified, energy independent) PSF.

    The effective area is essentially a linear spline, with built-in facilities
    for evaluation and plotting.

    Arguments
    ---------
    psf_file_path : str
        The path to the .psf FITS file containing the effective area table.

    rmax : float
        The maximum radial distance (in arcsec) for the PSF.
    """

    def __init__(self, psf_file_path, rmax=150.):
        """Constructor.
        """
        logger.info('Reading PSF data from %s...' % psf_file_path)
        hdu_list = fits.open(psf_file_path)
        hdu_list.info()
        _data = hdu_list['PSF'].data
        W = _data['W']
        sigma = _data['SIGMA']
        N = _data['N']
        r_c = _data['R_C']
        eta = _data['ETA']
        _x = numpy.linspace(0, rmax, 100)
        _y = W*numpy.exp(-(_x**2/(2*sigma**2))) + N*(1 + (_x/r_c)**2)**(-eta)
        fmt = dict(xname='r', xunits='arcsec', yname='PSF', yunits='sr$^-1$')
        xInterpolatedUnivariateSplineLinear.__init__(self, _x, _y, **fmt)

    def plot(self, num_points=1000, overlay=False, logx=False, logy=True,
             show=True):
        """Overloaded plot method (with default log scale on the y-axis).
        """
        xInterpolatedUnivariateSplineLinear.plot(self, num_points, overlay,
                                                 logx, logy, show)


def main():
    """
    """
    import os
    from ximpol import XIMPOL_IRF

    file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.psf')
    psf = xPointSpreadFunction(file_path)
    psf.plot()


if __name__ == '__main__':
    main()
