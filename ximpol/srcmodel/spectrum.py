#!/urs/bin/env python
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
from ximpol.core.rand import xUnivariateGenerator, xUnivariateAuxGenerator
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.core.spline import xInterpolatedBivariateSplineLinear


class xSourceSpectrum(xInterpolatedBivariateSplineLinear):

    """Class representing a count spectrum, i.e., the convolution of the
    source photon spectrum and the detector effective area.
    """

    def __init__(self, dNdE, E, t):
        """Constructor.
        """
        fmt = dict(xname='Time', xunits='s', yname='Energy',
                   yunits='keV',  zname='dN/dE',
                   zunits='cm$^{-2}$ s$^{-1}$ keV$^{-1}$')
        self.dNdE = dNdE
        xInterpolatedBivariateSplineLinear.__init__(self, t, E, dNdE, **fmt)
        self.light_curve = self.build_light_curve()

    def build_light_curve(self):
        """Build the light curve, i.e., a linear spline of the integral
        flux values as a function of time.
        """
        fmt = dict(rvname=self.xname, rvunits=self.xunits,
                   pdfname='Light curve', pdfunits='cm$^{-2}$ s$^{-1}$')
        _f = numpy.array([self.vslice(_t).norm() for _t in self.x])
        return xUnivariateGenerator(self.x, _f, **fmt)


def power_law(C, Gamma):
    """Photon energy spectrum as a function of energy and time.

    If C and Gamma are callable, we assume that the argument of the __call__
    function is the time, and this is how we treat them internally.
    """
    if hasattr(C, '__call__') and hasattr(Gamma, '__call__'):
        def _function(E, t):
            return C(t)*numpy.power(E, -Gamma(t))
    else:
        def _function(E, t):
            return C*numpy.power(E, -Gamma)
    return _function

def constant(C):
    """
    """
    def _function(E,t):
        return C
    return _function


class xTabulatedStationarySpectrum(xInterpolatedUnivariateSplineLinear):

    """Class describing a time-independent spectral model where the flux
    as a function of energy is tabulated at discrete values (e.g., read from
    a file).

    This is essentially a xInterpolatedUnivariateSplineLinear whose __call__
    method is overloaded to (i) provide a different signature, where one
    can pass two arguments, energy and time, although the second is not
    used; and (ii) allow to evaluate the spline on multidimensional vectors such
    as meshgrids, which is what happens when the spectrum is folded with
    the effective area to create the count spectrum.
    """

    def __init__(self, energy, flux, xname=None, xunits=None, yname=None,
                 yunits=None, optimize=False, tolerance=1e-4):
        """Constructor.
        """
        if xname is None:
            xname = 'Energy'
        if xunits is None:
            xunits='keV'
        if yname is None:
            yname='Flux'
        if yunits is None:
            yunits='cm$^{-2}$ s$^{-1}$'
        xInterpolatedUnivariateSplineLinear.__init__(self, energy, flux, xname,
                                                     xunits, yname, yunits,
                                                     optimize, tolerance)

    def __call__(self, E, t):
        """Overloaded __call__ method.

        Warning
        -------
        This is a major interface change with respect to the underlying scipy
        class. Also, there might be ways to be more efficient, if the only
        use case is to evaluate the spline over meshgrids, i.e., call the thing
        on a row and then tile the result.
        """
        if E.ndim > 1:
            y = xInterpolatedUnivariateSplineLinear.__call__(self, E.flatten())
            return numpy.reshape(y, E.shape)
        return xInterpolatedUnivariateSplineLinear.__call__(self, E)


class xCountSpectrum(xUnivariateAuxGenerator):

    """Class representing a count spectrum, i.e., the convolution of the
    source photon spectrum and the detector effective area.
    """

    def __init__(self, dNdE, aeff, t):
        """Constructor.
        """
        fmt = dict(auxname='Time', auxunits='s', rvname='Energy',
                   rvunits='keV',  pdfname='dN/dE $\\times$ aeff',
                   pdfunits='s$^{-1}$ keV$^{-1}$')

        def _pdf(E, t):
            """Return the convolution between the effective area and
            the input photon spectrum.
            """
            return dNdE(E, t)*numpy.tile(aeff(E[0]), (t.shape[0], 1))

        xUnivariateAuxGenerator.__init__(self, t, aeff.x, _pdf, **fmt)
        self.light_curve = self.build_light_curve()

    def build_light_curve(self):
        """Build the light curve, i.e., a linear spline of the integral
        flux values as a function of time.
        """
        fmt = dict(rvname=self.xname, rvunits=self.xunits,
                   pdfname='Light curve', pdfunits='s$^{-1}$')
        _f = numpy.array([self.slice(_t).norm() for _t in self.x])
        return xUnivariateGenerator(self.x, _f, **fmt)





def main():
    """
    """
    import os
    from ximpol import XIMPOL_IRF
    from ximpol.irf.arf import xEffectiveArea

    def dNdE(E, t):
        """Function defining a time-dependent energy spectrum.
        """
        return 10.0*(1.0 + numpy.cos(t))*numpy.power(E, (-2.0 + 0.01*t))

    file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.arf')
    aeff = xEffectiveArea(file_path)
    t = numpy.linspace(0, 25, 100)
    c = xCountSpectrum(dNdE, aeff, t)
    c.light_curve.plot()
    c.plot()


if __name__ == '__main__':
    main()
