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


"""Custom random generation classes, mostly based on splines.
"""


import numpy

from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.core.spline import xInterpolatedBivariateSplineLinear


class xUnivariateGenerator(xInterpolatedUnivariateSplineLinear):

    """Univariate random number generator based on a linear interpolated
    spline.

    Args
    ----
    rv : array
        Array of points sampling the values of the random variable.

    pdf : array
        pdf values at the array rv.

    rvname: str, optional
        The name of the random variable.

    rvunits: str, optional
        The units for the random variable.

    pdfname: str, optional
        The name of the pdf.

    pdfunits: str, optional
        The units for the pdf.
    """

    def __init__(self, rv, pdf, rvname=None, rvunits=None, pdfname='pdf',
                 pdfunits=None):
        """ Constructor.
        """
        if pdfunits is None and rvunits is not None:
            pdfunits = '1/%s' % rvunits
        xInterpolatedUnivariateSplineLinear.__init__(self, rv, pdf, rvname,
                                                     rvunits, pdfname, pdfunits)
        self.ppf = self.build_ppf()

    def pdf(self, rv):
        """Return the pdf value(s) at the point(s) rv.
        """
        return self(rv)

    def rvs(self, size=1):
        """Return random variates of arbitrary size.
        """
        return self.ppf(numpy.random.sample(size))


class xUnivariateAuxGenerator(xInterpolatedBivariateSplineLinear):

    """Univariate random generator where the pdf of the random variable
    might depend on an additional auxiliary variable.

    Internally, a meshgrid of (rv, aux) values is created, and the pdf
    is evaluated on the grid to construct an interpolated linear bivariate
    spline. A vertical slice of the bivariate spline at any given value of the
    auxiliary variable is the actual pdf of the random variable at that value of
    the auxiliary variable.

    Args
    ----
    aux : array
        Array of points sampling the values of the auxiliary variable (does\
        not need to be the same shape of rv).

    rv : array
        Array of points sampling the values of the random variable.

    pdf : callable or array
        A callable expressing the pdf as a function of a given pair of values\
        of rv and aux (in this order).

    auxname: str, optional
        The name of the auxiliary variable.

    auxunits: str, optional
        The units for the auxiliary variable.

    rvname: str, optional
        The name of the random variable.

    rvunits: str, optional
        The units for the random variable.

    pdfname: str, optional
        The name of the pdf.

    pdfunits: str, optional
        The units for the pdf.

    Example
    -------
    Suppose you have a power-law photon spectrum whose normalization and
    spectral index depend on time through the (bogus) relations

    >>> C(t) = 10*(1 + cos(t))
    >>> Gamma(t) = -2.0 + 0.01*t.

    We can construct a function encapsulating the spectrum by doing, say

    >>> def dNdE(E, t):
    >>>     \"""Function defining a time-dependent energy spectrum.
    >>>     \"""
    >>>     return 10.0*(1.0 + numpy.cos(t))*numpy.power(E, (-2.0 + 0.01*t))

    (note that we're using the numpy native functions, so that our callable
    can handle numpy arrays and evaluate the callable in an arbitrary number
    of points via broadcast). At this point we can sample the time (our
    auxiliary variable) and the energy (our random variable) on arbitrary grids
    of points and construct a random generator which is aware of the auxiliary
    variable:

    >>> from ximpol.core.rand import xUnivariateAuxGenerator
    >>>
    >>> t = numpy.linspace(0, 100, 100)
    >>> E = numpy.linspace(1, 10, 100)
    >>> fmt = dict(auxname='Time', auxunits='s', rvname='Energy', rvunits='keV',  pdfname='dN/dE')
    >>> gen = xUnivariateAuxGenerator(t, E, dNdE, **fmt)

    Given a vector of times, our generator will return a vector of
    (pseudo-random) energies, extracted with the spectral parameters that
    are appropriate for those times, e.g.:

    >>> t = numpy.random.uniform(0, 100, 1000000)
    >>> E = gen.rvs(t)

    Note
    ----
    `pdf` can be either a callable or an arraf shape (aux.size, rv.size).
    If `pdf` is a callable, than a meshgrid is created and the callable is
    evaluated on the meshgrid itself.
    """
    def __init__(self, aux, rv, pdf, auxname='aux', auxunits=None, rvname='rv',
                 rvunits=None, pdfname=None, pdfunits=None):
        """Constructor.
        """
        if pdfname is None:
            pdfname = 'pdf(%s; %s)' % (rvname, auxname)
        if pdfunits is None and rvunits is not None:
            pdfunits = '1/%s' % rvunits
        if hasattr(pdf, '__call__'):
            _rv, _aux = numpy.meshgrid(rv, aux)
            pdf = pdf(_rv, _aux)
        xInterpolatedBivariateSplineLinear.__init__(self, aux, rv, pdf,
                                                    auxname, auxunits,
                                                    rvname, rvunits,
                                                    pdfname, pdfunits)
        self.vppf = self.build_vppf()

    def pdf(self, aux, rv):
        """Return the pdf value(s) at the point(s) (rv, aux).
        """
        return self(aux, rv)

    def slice(self, aux):
        """Return the one-dimensional pdf for a given value of the auxiliary
        variable.
        """
        return self.vslice(aux)

    def rvs(self, aux):
        """Return random variates for a given array of values of the auxiliary
        variable.
        """
        return self.vppf(aux, numpy.random.sample(len(aux)))


def main():
    """
    """
    def dNdE(E, t):
        """Function defining a time-dependent energy spectrum.
        """
        return 10.0*(1.0 + numpy.cos(t))*numpy.power(E, (-2.0 + 0.01*t))

    t = numpy.linspace(0, 100, 100)
    E = numpy.linspace(1, 10, 100)
    fmt = dict(auxname='Time', auxunits='s', rvname='Energy', rvunits='keV',
               pdfname='dN/dE')
    gen = xUnivariateAuxGenerator(t, E, dNdE, **fmt)
    gen.plot()
    t = numpy.random.uniform(0, 100, 10)
    print (gen.rvs(t))


if __name__ == '__main__':
    main()
