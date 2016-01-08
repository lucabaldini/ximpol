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

    Example
    -------
    >>> from ximpol import XIMPOL_IRF
    >>> from ximpol.utils.matplotlib_ import pyplot as plt
    >>>
    >>> file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.psf')
    >>> psf = xPointSpreadFunction(file_path)
    >>> print(psf.rvs(10))
    >>> ra, dec = 5.0, 12.3
    >>> print(psf.smear_single(ra, dec, 10))

    Note
    ----
    The parametrization is taken from `Fabiani et al., 2014
    <http://arxiv.org/abs/1403.7200>`_, table 2. The PSF is technically
    energy dependent, but the dependence is not wild and for the moment
    we stick with the values at 4.51 keV in the table.
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
        self.ppf = self.build_ppf()

    def plot(self, num_points=1000, overlay=False, logx=False, logy=True,
             show=True):
        """Overloaded plot method (with default log scale on the y-axis).
        """
        xInterpolatedUnivariateSplineLinear.plot(self, num_points, overlay,
                                                 logx, logy, show)

    def rvs(self, size=1):
        """Random variates.
        """
        return self.ppf(numpy.random.sample(size))

    def delta(self, size=1):
        """Return an array of random offset (in ra, dec or L, B) due to the PSF.

        The output is converted in degrees.

        Warning
        -------
        Is this formula correct?
        """
        rho = self.rvs(size)/3600.
        phi = numpy.random.uniform(0, 2*numpy.pi, size)
        return  rho*numpy.cos(phi), rho*numpy.sin(phi)

    def smear_single(self, ra, dec, num_times=1):
        """Smear a pair of coordinates.
        """
        delta_ra, delta_dec = self.delta(size=num_times)
        return ra + delta_ra, dec + delta_dec

    def smear_array(self, ra, dec):
        """Smear a pair of arrays of coordinates.
        """
        assert(ra.size == dec.size)
        delta_ra, delta_dec = self.delta(ra.size())
        return ra + delta_ra, dec + delta_dec


def main():
    """
    """
    import os
    from ximpol import XIMPOL_IRF
    from ximpol.utils.matplotlib_ import pyplot as plt

    file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.psf')
    psf = xPointSpreadFunction(file_path)
    print(psf.rvs(10))
    ra, dec = 5.0, 12.3
    print(psf.smear_single(ra, dec, 10))
    xmax = 50
    plt.hist(psf.rvs(100000), xmax, (0, xmax), rwidth=1, histtype='step',
             lw=2, normed=True)
    _x = numpy.linspace(0, xmax, xmax)
    _y = psf(_x)/psf.norm()
    plt.plot(_x, _y)
    plt.show()


if __name__ == '__main__':
    main()
