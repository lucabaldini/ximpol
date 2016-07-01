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
from ximpol.irf.base import OGIP_HEADER_SPECS
from ximpol.core.fitsio import xBinTableHDUBase
from ximpol.core.spline import xInterpolatedUnivariateSpline
from ximpol.core.rand import xUnivariateGenerator


class xBinTableHDUPSF(xBinTableHDUBase):

    """Binary table for the PSF extension of a psf file.
    """

    NAME = 'PSF'
    DATA_SPECS = [
        ('W'    , 'E', '1/sr'),
        ('SIGMA', 'E', 'arcsec'),
        ('N'    , 'E', '1/sr'),
        ('R_C'  , 'E', 'arcsec'),
        ('ETA'  , 'E')
    ]


class xPointSpreadFunction(xInterpolatedUnivariateSpline):

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

    MAX_RADIUS = 150.

    def __init__(self, psf_file_path):
        """Constructor.
        """
        logger.info('Reading PSF data from %s...' % psf_file_path)
        self.hdu_list = fits.open(psf_file_path)
        self.hdu_list.info()
        _data = self.hdu_list['PSF'].data
        W = _data['W']
        sigma = _data['SIGMA']
        N = _data['N']
        r_c = _data['R_C']
        eta = _data['ETA']
        # Tabulate the actual PASF values.
        _r = numpy.linspace(0, self.MAX_RADIUS, 250)
        _y = W*numpy.exp(-(_r**2/(2*sigma**2))) + N*(1 + (_r/r_c)**2)**(-eta)
        fmt = dict(xname='r', xunits='arcsec', yname='PSF', yunits='sr$^{-1}$')
        xInterpolatedUnivariateSpline.__init__(self, _r, _y, k=2, **fmt)
        # And now include the solid angle for the actual underlying random
        # generator.
        _y *= 2*numpy.pi*_r
        fmt = dict(rvname='r', rvunits='arcsec',
                   pdfname='$2 \\pi r \\times$ PSF', pdfunits='')
        self.generator = xUnivariateGenerator(_r, _y, k=1, **fmt)
        self.eef = self.build_eef()
        #print 2*numpy.pi*W*sigma**2 + numpy.pi*r_c**2*N/(eta - 1.)

    def build_eef(self):
        """
        """
        return self.generator.build_cdf()

    def plot(self, show=True):
        """Overloaded plot method (with default log scale on the y-axis).
        """
        from ximpol.utils.matplotlib_ import pyplot as plt
        plt.figure('PSF')
        xInterpolatedUnivariateSpline.plot(self, logy=True, show=False)
        plt.figure('Solid-angle convolution')
        self.generator.plot(logy=True, show=False)
        plt.figure('EEF')
        self.eef.plot(show=False)
        if show:
            plt.show()

    def rvs(self, size):
        """
        """
        return self.generator.rvs(size)

    def delta(self, size=1):
        """Return an array of random offset (in ra, dec or L, B) due to the PSF.

        Note the output is converted in degrees.
        """
        rho = self.rvs(size)/3600.
        phi = numpy.random.uniform(0, 2*numpy.pi, size)
        return rho*numpy.cos(phi), rho*numpy.sin(phi)

    def smear_single(self, ra, dec, num_times=1):
        """Smear a pair of coordinates for an arbitrary number of times.
        """
        delta_ra, delta_dec = self.delta(size=num_times)
        return ra + delta_ra/numpy.cos(numpy.radians(dec)), dec + delta_dec

    def smear(self, ra, dec):
        """Smear a pair of arrays of coordinates.
        """
        assert(ra.size == dec.size)
        delta_ra, delta_dec = self.delta(ra.size)
        return ra + delta_ra/numpy.cos(numpy.radians(dec)), dec + delta_dec


def main():
    """
    """
    import os
    from ximpol import XIMPOL_IRF
    from ximpol.utils.matplotlib_ import pyplot as plt

    file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.psf')
    psf = xPointSpreadFunction(file_path)
    print(psf.rvs(10))
    ra, dec = 1., 1.
    print(psf.smear_single(ra, dec, 10))
    rmax = 100
    plt.hist(psf.rvs(1000000), rmax, (0, rmax), rwidth=1, histtype='step',
             lw=2, normed=True)
    _r = numpy.linspace(0, rmax, rmax)
    _psf = psf(_r)/psf.norm()
    plt.plot(_r, _psf)
    plt.show()


if __name__ == '__main__':
    main()
