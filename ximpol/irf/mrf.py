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


class xAzimuthalResponseGenerator(xUnivariateAuxGenerator):

    """Random number generator for the azimuthal response of the polarimeter.

    Here is the basic underlying math. Typically the response of a polarimeter
    to monochromatic, linearly polarized incident radiation is of the form:

    .. math:: N(\\phi) = A + B \\cos^2(\\phi - \\phi_0).

    This can be conveniently rewritten in term of the overall normalization
    (i.e., the total number of events) and the visibility of the modulation,
    defined as

    .. math::
        \\xi = \\frac{N_\\text{max} - N_\\text{min}}
        {N_\\text{max} + N_\\text{min}} = \\frac{B}{2A + B}

    (the visibility can also be characterized as the fraction of modulated
    events in a given distribution, as can be readily demonstrated.) The
    angular response now becomes

    .. math::
        N(\\phi) = \\frac{N_0}{\\pi} \\left[
        \\frac{(1 - \\xi)}{2} + \\xi \\cos^2(\\phi - \\phi_0)
        \\right].

    In terms of throwing random numbers, the phase is a trivial constant that
    can be added after the fact (modulo 2pi), so effectively the
    relevant probability density function is

    .. math::
        \\text{pdf}(\\phi) = \\frac{1}{\\pi} \\left[
        \\frac{(1 - \\xi)}{2} + \\xi \\cos^2(\\phi) \\right],

    (with the visibility being our auxiliary variable) and the corresponding
    cumulative distribution function is

    .. math::
        \\text{cdf}(\\phi) = \\frac{1}{2\\pi} \\left(
        \\phi + \\frac{\\xi}{2}\\sin{(2\\phi)} \\right).

    It is unfortunate that the cdf cannot be inverted, otherwise we would
    have no need to interpolate for generating random numbers according to
    this distribution.
    """

    def __init__(self):
        """Constructor.
        """
        _visibility = numpy.linspace(0., 1., 100)
        _phi = numpy.linspace(0., 2*numpy.pi, 100)
        fmt = dict(auxname='$\\xi$', rvname='$\\phi$', rvunits='rad')
        xUnivariateAuxGenerator.__init__(self, _visibility, _phi, self.pdf,
                                         **fmt)

    @classmethod
    def pdf(self, phi, visibility):
        """Evaluate the underlying one-dimensional pdf for a given value of the
        visibility, and assuming that the phase of the modulation is zero.

        Arguments
        ---------
        phi : float or array
            The (independent) azimuthal angle variable, in radians.

        visibility : float or array
            The visibility of the modulation, in [0--1].
        """
        return (0.5*(1. - visibility) +\
                visibility*numpy.power(numpy.cos(phi), 2.0))/numpy.pi

    @classmethod
    def cdf(self, phi, visibility):
        """Evaluate the underlying one-dimensional cdf for a given value of the
        visibility, and assuming that the phase of the modulation is zero.

        Warning
        -------
        We could overload the build_vpppf method for the class using this,
        since we have an analytic expression.
        """
        return (phi + 0.5*visibility*numpy.sin(2.*phi))/(2*numpy.pi)

    def rvs_phi(self, visibility, phase):
        """Generate random variates for any visibility and phase values.

        This is essentially calling the underlying xUnivariateAuxGenerator.rvs()
        method passing the visibility array as an argument and adding the phase
        array modulo 2pi.
        """
        return numpy.mod(self.rvs(visibility) + phase, 2*numpy.pi)

    @classmethod
    def fit_function(self, phi, visibility, phase, normalization):
        """Convenience function (with the phase back in) to allow histogram
        fitting.
        """
        return normalization*self.pdf((phi - phase), visibility)

    @classmethod
    def fit_histogram(self, histogram, fit_normalization=False):
        """Fit an azimuthal histogram.
        """
        from scipy.optimize import curve_fit
        _y, binning, patches = histogram
        _x = 0.5*(binning[1:] + binning[:-1])
        norm = _y.sum()*(binning[1] - binning[0])
        if not fit_normalization:
            p0 = (0.5, 0.5)
            # Wrap the fit function, keeping the normalization frozen.
            def f(phi, visibility, phase):
                return self.fit_function(phi, visibility, phase, norm)
        else:
            p0 = (0.5, 0.5, norm)
            f = self.fit_function
        popt, pcov = curve_fit(f, _x, _y, p0, numpy.sqrt(_y))
        if not fit_normalization:
            # Add back the normalization to the parameter vector and covariance
            # matrix.
            popt = numpy.concatenate((popt, numpy.array([norm])))
            pcov = numpy.vstack((pcov, numpy.array([[0., 0.]])))
            pcov = numpy.hstack((pcov, numpy.array([[0., 0., 0.]]).transpose()))
        _xs = numpy.linspace(0, 2*numpy.pi, 250)
        _ys = self.fit_function(_xs, *popt)/(binning[1] - binning[0])
        spline = xInterpolatedUnivariateSplineLinear(_xs, _ys)
        mask = _y > 0.
        obs = _y[mask]
        exp = []
        for _min, _max in zip(binning[:-1], binning[1:]):
            exp.append(spline.integral(_min, _max))
        exp = numpy.array(exp)
        chisquare = ((exp - obs)**2/exp).sum()
        return xModulationFitResults(popt, pcov, chisquare, len(mask))


class xModulationFitResults:

    """Small convenience class encapsulating the result of a fit to an
    azimuthal angle distribution.

    This includes facilities for plotting and annotating the best-fit
    model (e.g., overlaying it onto the underlying fitted histogram).
    """

    def __init__(self, popt, pcov, chisquare, ndof):
        """Constructor.
        """
        self.popt = popt
        self.pcov = pcov
        self.chisquare = chisquare
        self.ndof = ndof
        self.visibility, self.phase, self.normalization = popt
        self.visibility_error, self.phase_error,\
            self.normalization_error = numpy.sqrt(pcov.diagonal())

    def plot(self, show=False, stat=True, *options):
        """Plot the fit results.
        """
        from ximpol.utils.matplotlib_ import pyplot as plt
        _x = numpy.linspace(0., 2*numpy.pi, 100)
        _y = xAzimuthalResponseGenerator.fit_function(_x, *self.popt)
        plt.plot(_x, _y, *options)
        if stat:
            plt.text(0.05, 0.05, self.latex(), transform=plt.gca().transAxes)
        if show:
            plt.show()

    def latex(self):
        """LaTeX formatting.
        """
        return '$\\xi = %.3f \\pm %.3f$, $\\phi = (%.2f \\pm %.2f)^\\circ$,'\
            '$\\chi^2/{\\rm ndf} = %.1f/%d$' %\
            (self.visibility, self.visibility_error, numpy.degrees(self.phase),
             numpy.degrees(self.phase_error), self.chisquare, self.ndof)

    def __str__(self):
        """String formatting.
        """
        return 'Visibility = %.3f +/- %.3f, phase = %.2f +/- %.2f deg' %\
            (self.visibility, self.visibility_error, numpy.degrees(self.phase),
             numpy.degrees(self.phase_error))


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

        Warning
        -------
        Obsolete, to be removed.
        """
        return B/(2*A + B)

    @classmethod
    def A(self, mu):
        """Return the `A` parameter for the azimuthal cos-square distribution
        given the modulation factor `mu`.

        (This is calculated in a such a way the overall distribution is
        normalized to 1.)

        Warning
        -------
        Obsolete, to be removed.
        """
        return (1 - mu)/(2*numpy.pi)

    @classmethod
    def B(self, mu):
        """Return the `B` parameter for the azimuthal cos-square distribution
        given the modulation factor `mu`.

        (This is calculated in a such a way the overall distribution is
        normalized to 1.)

        Warning
        -------
        Obsolete, to be removed.
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

        Warning
        -------
        Obsolete, to be removed.
        """
        return norm*(self.A(visibility) + self.B(visibility)*numpy.power(
            numpy.cos(numpy.radians(phi - phi0)), 2.))

    @classmethod
    def fit_histogram(self, histogram):
        """Fit an azimuthal histogram.

        Warning
        -------
        Obsolete, to be removed.
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
