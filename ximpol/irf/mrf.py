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
from ximpol.core.rand import xUnivariateAuxGenerator
from ximpol.core.spline import optimize_grid_linear


"""Header specifications for the MODFRESP extension of .mrf FITS files.
"""
MODFRESP_HEADER_SPECS = [
    ('EXTNAME' , 'MODFRESP', 'name of this binary table extension'),
    ('HDUCLAS1', 'RESPONSE', 'dataset relates to spectral response'),
    ('HDUCLAS2', 'MODFRESP', 'dataset contains modulation response')
] + OGIP_HEADER_SPECS


class xColDefsMODFRESP(xColDefsBase):

    """ximpol.irf.base.xColDefsBase subclass for the MODFRESP extension
    of .mrf FITS files.
    """

    COLUMN_SPECS = [
        ('ENERG_LO', 'E', 'keV'),
        ('ENERG_HI', 'E', 'keV'),
        ('MODFRESP', 'E', None)
    ]



class xModulationFitResults:

    """Small convenience class encapsulating the result of a fit to an
    azimuthal angle distribution.
    """

    def __init__(self, popt, pcov):
        """Constructor.
        """
        self.popt = popt
        self.pcov = pcov
        self.vis, self.phi0, self.norm = popt
        self.vis_err, self.phi0_err, self.norm_err = numpy.sqrt(pcov.diagonal())

    def plot(self, *options):
        """Plot the fit results.
        """
        from ximpol.utils.matplotlib_ import pyplot as plt
        x = numpy.linspace(0., 360., 100)
        y = xModulationFactor.modulation_function(x, *self.popt)
        plt.plot(x, y, *options)

    def __str__(self):
        """String formatting.
        """
        return 'Visibility = %.3f +/- %.3f, angle = %.2f +/- %.2f' %\
            (self.vis, self.vis_err, self.phi0, self.phi0_err)


class xModulationFactor(xInterpolatedUnivariateSplineLinear):

    """Class describing the modulation factor.

    The effective area is essentially a linear spline, with built-in facilities
    for evaluation and plotting.

    Arguments
    ---------
    mrf_file_path : str
        The path to the .mrf FITS file containing the effective area table.

    Example
    -------
    >>> import os
    >>> import numpy
    >>> from ximpol import XIMPOL_IRF
    >>>
    >>> file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.mrf')
    >>> modf = xModulationFactor(file_path)
    >>> x = numpy.arange(1, 10, 1)
    >>> print(modf(x))
    >>> modf.plot(overlay=False)
    """

    def __init__(self, mrf_file_path):
        """Constructor.
        """
        logger.info('Reading modulation factor data from %s...' % mrf_file_path)
        hdu_list = fits.open(mrf_file_path)
        hdu_list.info()
        _data = hdu_list['MODFRESP'].data
        _x = 0.5*(_data.field('ENERG_LO') + _data.field('ENERG_HI'))
        _y = _data.field('MODFRESP')
        hdu_list.close()
        fmt = dict(xname='Energy', xunits='keV', yname='Modulation factor',
                   optimize=True, tolerance=1e-4)
        xInterpolatedUnivariateSplineLinear.__init__(self, _x, _y, **fmt)

    @classmethod
    def mu(self, A, B):
        """Return the modulation factor for an azimuthal distribution like

        >>> N(phi) = A + B*(cos(phi - phi0))**2

        with `A` and `B` passed as arguments.
        """
        return B/(2*A + B)

    @classmethod
    def A(self, mu):
        """Return the `A` parameter for the azimuthal cos-square distribution
        given the modulation factor `mu`.

        (This is calculated in a such a way the overall distribution is
        normalized to 1.)
        """
        return (1 - mu)/(2*numpy.pi)

    @classmethod
    def B(self, mu):
        """Return the `B` parameter for the azimuthal cos-square distribution
        given the modulation factor `mu`.

        (This is calculated in a such a way the overall distribution is
        normalized to 1.)
        """
        return mu/numpy.pi

    @classmethod
    def modulation_function(self, phi, visibility, phi0, norm=1):
        """The internal representation of the azimuthal cos-square distribution.

        This is used to construct the underlying generator, but can be used
        to fit modulation histograms, as well.

        Arguments
        ---------
        phi : float or array
            The value(s) of the azimuthal angle where we want to evaluate the\
            modulation function.

        mu : float
            The visibility of the modulation, i.e., (max - min)/(max + min).

        angle : float
            The modulation angle in degrees [0--360].
        """
        return norm*(self.A(visibility) + self.B(visibility)*numpy.power(
            numpy.cos(numpy.radians(phi - phi0)), 2.))

    @classmethod
    def fit_histogram(self, histogram):
        """Fit an azimuthal histogram.
        """
        from scipy.optimize import curve_fit
        _y, binning, patches = histogram
        _x = 0.5*(binning[1:] + binning[:-1])
        pstart = (0.5, 0., _y.sum())
        popt, pcov = curve_fit(self.modulation_function, _x, _y, pstart,
                               numpy.sqrt(_y))
        return xModulationFitResults(popt, pcov)

    def build_generator(self, polarization_angle=0., polarization_degree=1.):
        """Construct the underlying generator to throw random numbers
        according to the proper distribution.

        Arguments
        ---------
        polarization_angle : float
            The polarization angle in degrees [0--360].

        polarization_degree : float
            The polarization degree [0--1].

        Warning
        -------
        More work is needed to fully support energy- and time-dependent
        polarization degree and agle. (We probably need some kind of
        trivariate spline data sctructure.)
        """
        _x = self.x.copy()
        _y = numpy.linspace(0., 360., 100)
        # Maybe the loop for _z could be done in a more numpy-like fashion?
        _z = numpy.zeros(shape = (_x.size, _y.size))
        for i, _xp in enumerate(_x):
            mu = self(_xp)
            for j, _yp in enumerate(_y):
                _z[i, j] = self.modulation_function(_yp, mu*polarization_degree,
                                                    polarization_angle)
        fmt = dict(auxname='Energy', auxunits='keV', rvname='Azimuthal angle',
                   rvunits='$^\circ$')
        self.generator = xUnivariateAuxGenerator(_x, _y, _z, **fmt)

    def rvs(self, E):
        """Return random variates for a given array of values of energies.
        """
        return self.generator.rvs(E)


def main():
    """
    """
    import os
    import numpy
    from ximpol import XIMPOL_IRF

    file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.mrf')
    modf = xModulationFactor(file_path)
    x = numpy.linspace(1, 10, 10)
    print(modf(x))
    modf.plot(overlay=True)
    modf.build_generator(polarization_angle=20.)
    modf.generator.plot()
    modf.generator.slice(5).plot()


if __name__ == '__main__':
    main()
