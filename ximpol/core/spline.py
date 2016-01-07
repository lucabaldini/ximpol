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


"""Spline utility module, building on top of the scipy.interpolate modules.
"""


import numpy
from scipy.interpolate import UnivariateSpline, InterpolatedUnivariateSpline
from scipy.interpolate import RectBivariateSpline

from ximpol.utils.logging_ import logger


class xUnivariateSplineBase:

    """Base class for all the univariate spline classes.

    The basic idea is to keep track of the original arrays passed to the
    interpolator and to support arithmetic operations. We also allow the user
    to supply optional arguments to control the ranges and specify names and
    units for the quantities involved.

    Args
    ----
    x : array
        Input x values (assumed to be sorted).

    y : array
        Input y values.

    xname: str, optional
        The name of the quantity on the x-axis.

    xunits: str, optional
        The units for the x-axis.

    yname: str, optional
        The name of the quantity on the y-axis.

    yunits: str, optional
        The units for the y-axis.

    Note
    ----
    This is a do-nothing class to be subclassed and not instantiated
    directly.
    """

    def __init__(self, x, y, xname=None, xunits=None, yname=None, yunits=None):
        """Constructor.
        """
        assert(len(x) == len(y))
        self.x = x.copy()
        self.y = y.copy()
        self.xname = xname
        self.xunits = xunits
        self.yname = yname
        self.yunits = yunits

    def xmin(self):
        """Return the minimum of the underlying x-array.
        """
        return self.x[0]

    def xmax(self):
        """Return the maximum of the underlying x-array.
        """
        return self.x[-1]

    def __mul__(self, other):
        """Overloaded multiplication operator.
        """
        assert(self.__class__.__name__ == other.__class__.__name__)
        _x = numpy.union1d(self.x, other.x)
        _y = self(_x)*other(_x)
        return self.__class__(_x, _y)

    def __add__(self, other):
        """Overloaded sum operator.
        """
        assert(self.__class__.__name__ == other.__class__.__name__)
        _x = numpy.union1d(self.x, other.x)
        _y = self(_x) + other(_x)
        return self.__class__(_x, _y)

    def __len__(self):
        """Return the lenght of the arrays used to construct the spline.
        """
        return len(self.x)

    @classmethod
    def label(self, name, units=None):
        """Compose an axis label given a name and some units.
        """
        if units is None:
            return name
        else:
            return '%s [%s]' % (name, units)

    def xlabel(self):
        """Return the x-label for a plot.
        """
        return self.label(self.xname, self.xunits)

    def ylabel(self):
        """Return the y-label for a plot.
        """
        return self.label(self.yname, self.yunits)

    def norm(self):
        """Return the integral over the entire spline domain.
        """
        return self.integral(self.xmin(), self.xmax())

    def build_cdf(self):
        """Create the cumulative distribution function.
        """
        _x = self.x.copy()
        _y = numpy.array([self.integral(0, _xp) for _xp in _x])/self.norm()
        return self.__class__(_x, _y)

    def build_ppf(self):
        """Create the percent point function (or inverse of cdf).
        """
        _y = self.x.copy()
        _x = numpy.array([self.integral(0, _xp) for _xp in _y])
        _x, _mask = numpy.unique(_x, return_index=True)
        _x/= self.norm()
        _y = _y[_mask]
        return self.__class__(_x, _y)

    def plot(self, num_points=1000, overlay=False, logx=False, logy=False,
             show=True):
        """Plot the spline.

        Args
        ----
        num_points : int, optional
            The number of sampling points to be used to draw the spline.

        overlay : bool, optional
            If True, the original arrays passed to the spline are overlaid.

        show : bool, optional
            If True, `plt.show()` is called at the end, interrupting the flow.
        """
        from ximpol.utils.matplotlib_ import pyplot as plt
        _x = numpy.linspace(self.xmin(), self.xmax(), num_points)
        _y = self(_x)
        if overlay:
            plt.plot(_x, _y, '-', self.x, self.y, 'o')
        else:
            plt.plot(_x, _y, '-')
        if self.xname is not None:
            plt.xlabel(self.xlabel())
        if self.yname is not None:
            plt.ylabel(self.ylabel())
        if logx:
            plt.gca().set_xscale('log')
        if logy:
            plt.gca().set_yscale('log')
        if show:
            plt.show()


class xUnivariateSpline(xUnivariateSplineBase, UnivariateSpline):

    """Light-weight wrapper over the scipy `UnivariateSpline
    <http://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html>`_
    class (see the corresponding documentation for the meaning of the
    parameters passed to the constructor).

    Note
    ----
    Note that the interface to the base class has changed from numpy 0.14.
    An `ext` argument can be passed to the constructor starting with scipy
    0.15 to control the extrapolation behavior and a `check_finite` argument is
    available in 0.16 to avoid `nans` in the input data.
    We currently do not use either one.
    """

    def __init__(self, x, y, w=None, bbox=[None, None], k=3, s=None,
                 xname=None, xunits=None, yname=None, yunits=None):
        """Constructor.
        """
        xUnivariateSplineBase.__init__(self, x, y, xname, xunits, yname, yunits)
        UnivariateSpline.__init__(self, self.x, self.y, w, bbox, k, s)


class xInterpolatedUnivariateSpline(xUnivariateSplineBase,
                                    InterpolatedUnivariateSpline):

    """Light-weight wrapper over the scipy `InterpolatedUnivariateSpline
    <http://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.InterpolatedUnivariateSpline.html>`_
    class (see the corresponding documentation for the meaning of the
    parameters passed to the constructor).

    Note
    ----
    Note that the interface to the base class has changed from numpy 0.14.
    An `ext` argument can be passed to the constructor starting with scipy
    0.15 to control the extrapolation behavior and a `check_finite` argument is
    available in 0.16 to avoid `nans` in the input data.
    We currently do not use either one.
    """

    def __init__(self, x, y, w=None, bbox=[None, None], k=3,
                 xname=None, xunits=None, yname=None, yunits=None):
        """Constructor.
        """
        xUnivariateSplineBase.__init__(self, x, y, xname, xunits, yname, yunits)
        InterpolatedUnivariateSpline.__init__(self, self.x, self.y, w, bbox, k)


class xInterpolatedUnivariateSplineLinear(xInterpolatedUnivariateSpline):

    """xInterpolatedUnivariateSplineLinear subclass implementing the simplest
    possible linear interpolator.

    Example
    -------
    >>> from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
    >>>
    >>> x = numpy.linspace(0, 2*numpy.pi, 100)
    >>> y = numpy.sin(x)
    >>> s = xInterpolatedUnivariateSplineLinear(x, y, xname='x', xunits='au')
    >>> s.plot()
    """

    def __init__(self, x, y, xname=None, xunits=None, yname=None, yunits=None):
        """ Constructor.
        """
        xInterpolatedUnivariateSpline.__init__(self, x, y, None, [None, None],
                                               1, xname, xunits, yname, yunits)


class xBivariateSplineBase:

    """Base class for all the bivariate spline classes.

    This is somewhat similar in spirit to the corresponding univariate base
    class, except that the additional functionalities are, for the moment,
    limited to bookkeeping and plotting facilities.

    Args
    ----
    x : array
        Input x values (assumed to be sorted).

    y : array
        Input y values (assumed to be sorted).

    z : array
        Input z values, with shape (x.size,y.size).

    xname: str, optional
        The name of the quantity on the x-axis.

    xunits: str, optional
        The units for the x-axis.

    yname: str, optional
        The name of the quantity on the y-axis.

    yunits: str, optional
        The units for the y-axis.

    zname: str, optional
        The name of the quantity on the z-axis.

    zunits: str, optional
        The units for the z-axis.

    Note
    ----
    This is a do-nothing class to be subclassed and not instantiated
    directly.
    """

    def __init__(self, x, y, z, xname=None, xunits=None, yname=None,
                 yunits=None, zname=None, zunits=None):
        """Constructor.
        """
        self.x = x
        self.y = y
        self.z = z
        self.xname = xname
        self.xunits = xunits
        self.yname = yname
        self.yunits = yunits
        self.zname = zname
        self.zunits = zunits

    def xmin(self):
        """Return the minimum of the underlying x-array.
        """
        return self.x[0]

    def xmax(self):
        """Return the maximum of the underlying x-array.
        """
        return self.x[-1]

    def ymin(self):
        """Return the minimum of the underlying y-array.
        """
        return self.y[0]

    def ymax(self):
        """Return the maximum of the underlying y-array.
        """
        return self.y[-1]

    def xlabel(self):
        """Return the x-label for a plot.
        """
        return xUnivariateSplineBase.label(self.xname, self.xunits)

    def ylabel(self):
        """Return the y-label for a plot.
        """
        return xUnivariateSplineBase.label(self.yname, self.yunits)

    def zlabel(self):
        """Return the z-label for a plot.
        """
        return xUnivariateSplineBase.label(self.zname, self.zunits)


class xInterpolatedBivariateSplineLinear(xBivariateSplineBase,
                                         RectBivariateSpline):

    """Bivariate linear interpolated spline on a rectangular grid.
    """

    def __init__(self, x, y, z, xname=None, xunits=None, yname=None,
                 yunits=None, zname=None, zunits=None):
        """Constructor.
        """
        xBivariateSplineBase.__init__(self, x, y, z, xname, xunits, yname,
                                      yunits, zname, zunits)
        RectBivariateSpline.__init__(self, x, y, z,
                                     bbox=[None, None, None, None],
                                     kx=1, ky=1, s=0)

    def __call__(self, x, y, dx=0, dy=0, grid=False):
        """Overloaded __call__method.

        Here we basically override the default value of the `grid` parameter
        from `True` to `False`, since we're typically interested in evaluating
        the splined at given physical coordinates, rather than grid points.
        """
        return RectBivariateSpline.__call__(self, x, y, None, dx, dy, grid)

    def vslice(self, x, num_points=1000):
        """Return a vertical slice at a given x of the bivariate spline.

        Args
        ----
        x : float
            The x value at which the vertical slice should be calculated.

        num_points : int, optional
            The number of sampling points for the output univariate spline.
        """
        _x = numpy.linspace(self.ymin(), self.ymax(), num_points)
        _y = self(x, _x)
        fmt = dict(xname=self.yname, xunits=self.yunits, yname=self.zname,
                   yunits=self.zunits)
        return xInterpolatedUnivariateSplineLinear(_x, _y, **fmt)

    def hslice(self, y, num_points=1000):
        """Return an horizontal slice at a given y of the bivariate spline.

        Args
        ----
        y : float
            The y value at which the horizontal slice should be calculated.

        num_points : int, optional
            The number of sampling points for the output univariate spline.
        """
        _x = numpy.linspace(self.xmin(), self.xmax(), num_points)
        _y = self(_x, y)
        fmt = dict(xname=self.xname, xunits=self.xunits, yname=self.zname,
                   yunits=self.zunits)
        return xInterpolatedUnivariateSplineLinear(_x, _y, **fmt)

    def plot(self, num_pointsx=100, num_pointsy=100, num_contours=75,
             show=True):
        """Plot the spline.

        Args
        ----
        num_pointsx : int
            The number of x sampling points to be used to draw the spline.

        num_pointsy : int
            The number of y sampling points to be used to draw the spline.

        num_contours : int
            The number of contours for the color plot.

        show : bool, optional
            If True, `plt.show()` is called at the end, interrupting the flow.
        """
        from ximpol.utils.matplotlib_ import pyplot as plt
        _x = numpy.linspace(self.xmin(), self.xmax(), num_pointsx)
        _y = numpy.linspace(self.ymin(), self.ymax(), num_pointsy)
        _x, _y = numpy.meshgrid(_x, _y)
        _z = self(_x, _y, grid=False)
        contour = plt.contourf(_x, _y, _z, num_contours)
        bar = plt.colorbar()
        if self.xname is not None:
            plt.xlabel(self.xlabel())
        if self.yname is not None:
            plt.ylabel(self.ylabel())
        if self.zname is not None:
            bar.set_label(self.zlabel())
        if show:
            plt.show()


def main():
    x = numpy.linspace(0, 2*numpy.pi, 100)
    y = numpy.sin(x)
    s = xInterpolatedUnivariateSplineLinear(x, y, 'x', 'au', 'y')
    s.plot()


if __name__ == '__main__':
    main()
