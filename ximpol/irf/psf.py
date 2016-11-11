#!/usr/bin/env python
#
# Copyright (C) 2015, the ximpol team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
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


def gauss_king(r, W, sigma, N, r_c, eta):
    """Functional representation of the Gaussian plus King PSF profile
    described in `Fabiani et al., 2014 <http://arxiv.org/abs/1403.7200>`_,
    equation (2):

    .. math::
        \\text{PSF}(r) = W \\exp^{-(\\frac{r^2}{2\\sigma^2})} +
        N\\left( 1 + \\left( \\frac{r}{r_c} \\right)^2 \\right)^{-\\eta}

    Arguments
    ---------
    r : float or array
        The radial distance from the true source position is arcsec.

    W : float
        Normalization of the Gaussian component.

    sigma : float
        Width of the Gaussian component.

    N : float
        Normalization of the King component.

    r_c : float
        Characteristic radius of the King component.

    eta : float
        Exponent of the King component.
    """
    return W*numpy.exp(-(r**2/(2*sigma**2))) + N*(1 + (r/r_c)**2)**(-eta)


def gauss_king_eef_at_infinity(W, sigma, N, r_c, eta):
    """Return the value of the Encircled Energy Fraction (EEF) at infinity,
    given the parameters of the functional representation, see equation (4)
    of `Fabiani et al., 2014 <http://arxiv.org/abs/1403.7200>`_.

    .. math::
        \\text{EEF}(\\infty) = 2\\pi W\\sigma^2 +
        \\pi\\frac{r_c^2 N}{\\eta - 1}

    Arguments
    ---------
    r : float or array
        The radial distance from the true source position is arcsec.

    W : float
        Normalization of the Gaussian component.

    sigma : float
        Width of the Gaussian component.

    N : float
        Normalization of the King component.

    r_c : float
        Characteristic radius of the King component.

    eta : float
        Exponent of the King component.
    """
    return 2*numpy.pi*W*sigma**2 + numpy.pi*r_c**2*N/(eta - 1.)


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
    PARAM_NAMES = ['W', 'sigma', 'N', 'r_c', 'eta']

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
        self.__params = (W, sigma, N, r_c, eta)
        # Tabulate the actual PSF values.
        _r = numpy.linspace(0, self.MAX_RADIUS, 250)
        _y = gauss_king(_r, *self.__params)
        fmt = dict(xname='r', xunits='arcsec', yname='PSF', yunits='sr$^{-1}$')
        xInterpolatedUnivariateSpline.__init__(self, _r, _y, k=2, **fmt)
        # Include the solid angle for the actual underlying random generator.
        _y *= 2*numpy.pi*_r
        fmt = dict(rvname='r', rvunits='arcsec',
                   pdfname='$2 \\pi r \\times$ PSF', pdfunits='')
        self.generator = xUnivariateGenerator(_r, _y, k=1, **fmt)
        # Finally, calculate the
        self.eef, self.hew = self.build_eef()
        logger.info(self)

    def build_eef(self):
        """Build the Encircled Energy Fraction (EEF) as a function of r.

        And, while we're at it, we also calculate and cache the HEW.
        """
        _r = self.x
        _y = numpy.array([self.generator.integral(_r[0], _rp) for _rp in _r])
        _y /= gauss_king_eef_at_infinity(*self.__params)
        hew = 2*xInterpolatedUnivariateSpline(_y, _r, k=2)(0.5)
        fmt = dict(xname='r', xunits='arcsec', yname='EEF')
        return xInterpolatedUnivariateSpline(_r, _y, k=1, **fmt), hew

    def view(self, show=True):
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
        """Extract values of the radial distance according to the PSF shape.
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

    def __str__(self):
        """String formatting.
        """
        text = 'Gauss + King PSF, '
        for i, name in enumerate(self.PARAM_NAMES):
            text += '%s = %.3e, ' % (name, self.__params[i])
        text += 'HEW = %.1f arcsec' % self.hew
        return text

    def draw_psf_circle(self, image, x, y, text='PSF', color='white', lw=2,
                        number=False):
        """Add the PSF circle to the image with labels. This function must be
        called after the (possible) image recenter.

        Note the x and y are coordinates relative to the figure axes (0.0 is
        left or bottom and 1.0 is right or top).
        """
        psf_rad = self.hew/(2.*60.**2) #degrees
        xpix_low, xpix_high = image._ax1.get_xbound()
        ypix_low, ypix_high = image._ax1.get_ybound()
        xpix_circ = x*(xpix_high-xpix_low) + xpix_low
        ypix_circ = y*(ypix_high-ypix_low) + ypix_low
        ra_circ, dec_circ = image.pixel2world(xpix_circ, ypix_circ)
        dec_text_up = dec_circ+1.7*psf_rad
        dec_text_down = dec_circ-1.9*psf_rad
        image.add_label(ra_circ, dec_text_up, text, size='x-large', color=color,
                        horizontalalignment='center')
        image.show_circles(ra_circ, dec_circ, psf_rad, lw=lw, color=color)
        if number:
            text_psf = '%d"' %round(self.hew)
            image.add_label(ra_circ, dec_text_down, text_psf, size='large',
                           color=color, horizontalalignment='center')
