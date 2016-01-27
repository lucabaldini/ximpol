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
    events in a given distribution, as can be readily demonstrated, and it is
    the quantity that the detector is effectively measuring).
    The angular response now becomes

    .. math::
        N(\\phi) = \\frac{N_0}{\\pi} \\left[
        \\frac{(1 - \\xi)}{2} + \\xi \\cos^2(\\phi - \\phi_0)
        \\right].

    For completeness, the visibility, the modulation factor and the polarization
    degree for a monocromatic source are related to each other through:

    .. math::
        \\xi(E, t) = P(E, t) \\times \\mu(E)

    (i.e., at any given energy the modulation factor is the visibility of the
    modulation of the detector response for 100% linearly polarized incident
    radiation).

    In terms of throwing random numbers, the phase is a trivial constant that
    can be added after the fact (modulo 2pi), so effectively the
    relevant probability density function is

    .. math::
        \\text{pdf}(\\phi) = \\frac{1}{\\pi} \\left[
        \\frac{(1 - \\xi)}{2} + \\xi \\cos^2(\\phi) \\right],

    .. image:: ../figures/test_azimuthal_resp_pdf.png

    The corresponding cumulative distribution function is

    .. math::
        \\text{cdf}(\\phi) = \\frac{1}{2\\pi} \\left(
        \\phi + \\frac{\\xi}{2}\\sin{(2\\phi)} \\right),

    and it is unfortunate that it cannot be inverted, otherwise we would
    have no need to interpolate for generating random numbers according to
    this distribution.

    .. image:: ../figures/test_azimuthal_resp_cdf.png

    From the prospecive of the code, this generator is a standard
    `xUnivariateAuxGenerator` where the azimuthal angle is our
    random variable and the visbility is our auxiliary variable. For any given
    value of the visibility, a vertical slice is providing the corresponding
    one-dimensional pdf.

    .. image:: ../figures/test_azimuthal_resp_generator.png

    The class also provide facilities to fit a histogram to recover the
    underlying modulation visibility and phase.

    Example
    -------
    >>> import numpy
    >>> from ximpol.utils.matplotlib_ import pyplot as plt
    >>> from ximpol.irf.mrf import xAzimuthalResponseGenerator
    >>>
    >>> generator = xAzimuthalResponseGenerator()
    >>> visibility = numpy.zeros(1000000)
    >>> visibility.fill(0.5)
    >>> phase = numpy.radians(45.)
    >>> phi = generator.rvs_phi(visibility, phase)
    >>> hist = plt.hist(phi, bins=numpy.linspace(0, 2*numpy.pi, 100),
                        histtype='step', label='Random angles')
    >>> fit_results = generator.fit_histogram(hist)
    >>> fit_results.plot()
    >>> plt.show()

    .. image:: ../figures/test_azimuthal_resp_rvs.png
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

        Arguments
        ---------
        phi : float or array
            The (independent) azimuthal angle variable, in radians.

        visibility : float or array
            The visibility of the modulation, in [0--1].

        Warning
        -------
        We could overload the build_vpppf method for the class using this,
        since we have an analytic expression. (This function is effectively
        not used at the moment.)
        """
        return (phi + 0.5*visibility*numpy.sin(2.*phi))/(2*numpy.pi)

    def rvs_phi(self, visibility, phase):
        """Generate random variates for any visibility and phase values.

        This is essentially calling the underlying xUnivariateAuxGenerator.rvs()
        method passing the visibility array as an argument and adding the phase
        array (modulo 2pi) to the output.

        Arguments
        ---------
        visibility : array
            An array of visibility values. (The function returns an equal-length
            array of phi values.)

        phase : float or array
            The phase of the modulation. (This can either be a vector or an
            array of the same length as `visibility`.)
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
        # Horrible hack.
        if popt[0] < 0.:
            popt[0] = -popt[0]
            popt[1] += 0.5*numpy.pi
        if popt[1] < 0.:
            popt[1] += numpy.pi
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

    def plot(self, show=False, stat=True, label=None, *options):
        """Plot the fit results.
        """
        from ximpol.utils.matplotlib_ import pyplot as plt
        _x = numpy.linspace(0., 2*numpy.pi, 100)
        _y = xAzimuthalResponseGenerator.fit_function(_x, *self.popt)
        plt.plot(_x, _y, *options)
        if label is not None:
            plt.text(0.05, 0.1, label, transform=plt.gca().transAxes)
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
        The path to the .mrf FITS file containing the modulation response table.


    To zero-th order, an `xModulationFactor` instance is an object capable of
    evaluating itself in a point or in a grid of points, and capable of
    plotting itself.

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
    >>> modf.plot()

    More interestingly, it can generate random `phi` values, given a vector
    of event energies and corresponding vectors (or simple floats) representing
    the polarization degree and angle corresponding to the energies themselves.
    Internally, any `xModulationFactor` has an `xAzimuthalResponseGenerator`
    object and when the `xModulationFactor.rvs_phi()` method is called,
    the polarization degree is multiplied by the modulation factor of the
    detector, evaluated at the right energy, and converted into a visibility
    value, after which the underlying `xAzimuthalResponseGenerator.rvs_phi()`
    is called.

    Example
    -------
    >>> import os
    >>> import numpy
    >>> from ximpol import XIMPOL_IRF
    >>>
    >>> file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.mrf')
    >>> modf = xModulationFactor(file_path)
    >>> # Throw 100000 random energy values.
    >>> energy = numpy.random.uniform(1, 10, 100000)
    >>> # This will create an array of 100000 phi values where the visibility
    >>> # of the modulation is tracking the modulation factor of the polarimeter
    >>> # (the phase is constant at 45 degrees).
    >>> phi = modf.rvs_phi(energy, 1., numpy.radians(45))
    """

    def __init__(self, mrf_file_path):
        """Constructor.
        """
        logger.info('Reading modulation factor data from %s...' % mrf_file_path)
        self.hdu_list = fits.open(mrf_file_path)
        self.hdu_list.info()
        _data = self.hdu_list['MODFRESP'].data
        _x = 0.5*(_data.field('ENERG_LO') + _data.field('ENERG_HI'))
        _y = _data.field('MODFRESP')
        fmt = dict(xname='Energy', xunits='keV', yname='Modulation factor',
                   optimize=True, tolerance=1e-4)
        xInterpolatedUnivariateSplineLinear.__init__(self, _x, _y, **fmt)
        self.generator = xAzimuthalResponseGenerator()

    def rvs_phi(self, energy, polarization_degree, polarization_angle):
        """Return random variates for a given array of values of energy,
        polarization degree and polarization angle.

        Arguments
        ---------
        energy : array
            An array of energy values. (The function returns an equal-length
            array of phi values.)

        polarization_degree : array or float
            The polarization degree, in [0--1]. (This can either be a vector
            or an array of the same length as `energy`.)

        polarization_angle : array or float
            The polarization angle, in radians. (This can either be a vector or
            an array of the same length as `energy`.)

        """
        visibility = self(energy)*polarization_degree
        return self.generator.rvs_phi(visibility, polarization_angle)


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
    modf.generator.plot()
    modf.generator.slice(0.5).plot()


if __name__ == '__main__':
    main()
