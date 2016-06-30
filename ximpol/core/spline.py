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

from ximpol.utils.logging_ import logger, abort


def interpolate(xa, ya, xb, yb, x):
    """Simple two-point linear interpolation/extrapolation.
    """
    return ya + (yb - ya)/(xb - xa)*(x - xa)

def optimize_grid_linear(x, y, tolerance=1e-4):
    """Optimize a pair of (x, y) arrays for the corresponding spline
    definition.

    This loops over the input arrays and removes unnecessary data points
    to minimize the length of the arrays necessary to the spline definition.

    Args
    ----
    x : array
        The input x-array.

    y : array
        The input y-array.

    tolerance : float
        The maximum relative difference between the generic yi value and the\
        estrapolation of the two previous optimized data points for the point\
        i to be removed.
    """
    assert(len(x) == len(y))
    logger.info('Optimizing grid with %d starting points...' % len(x))
    # Start a new series with the first two points of the input arrays.
    _x = [x[0], x[1]]
    _y = [y[0], y[1]]
    # Loop over the points 3 ... (N - 1).
    for i, (_xi, _yi) in enumerate(zip(x, y)[2:-1]):
        # Extrapolate the last two points of the new series to xi and
        # see how far we are from the actual yi.
        delta = interpolate(_x[-2], _y[-2], _x[-1], _y[-1], _xi) - _yi
        if abs(delta/_yi) > tolerance:
            # If the difference is larger than the tolerance, add a point.
            # (This has the drawback that we tend to add pairs of point at
            # each change of slope.)
            _x.append(_xi)
            _y.append(_yi)
            # Interpolate the points last and (last - 2) to (last - 1).
            delta = interpolate(_x[-3], _y[-3], _x[-1], _y[-1], _x[-2]) - _y[-2]
            if abs(delta/_y[-2]) < tolerance:
                # If the penultimate point was not necessary, remove it.
                _x.remove(_x[-2])
                _y.remove(_y[-2])
    # Append the last point of the original array to the list.
    _x.append(x[-1])
    _y.append(y[-1])
    _x, _y = numpy.array(_x), numpy.array(_y)
    logger.info('Done, %d points remaining.' % len(_x))
    return _x, _y


class xUnivariateSplineBase:

    """Base class for all the univariate spline classes.

    The basic idea is to keep track of the original arrays passed to the
    interpolator and to support arithmetic operations. We also allow the user
    to supply optional arguments to control the ranges and specify names and
    units for the quantities involved.

    Args
    ----
    x : array
        Input x values.

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
        # Make sure the input vectors have the same lengths.
        assert len(x) == len(y)
        # numpy.unique is returning a sorted copy of the unique values of the
        # x arrays, so this is effectively sorting x.
        self.x, _index = numpy.unique(x, return_index=True)
        # If some of the values were not unique, give up.
        assert len(self.x) == len(x)
        # Need to grab the y in the right order.
        self.y = y[_index]
        self.xname = xname
        self.xunits = xunits
        self.yname = yname
        self.yunits = yunits

    def dist(self, x, y):
        """
        """
        return abs((self(x) - y)/y)

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

    def __div__(self, other):
        """Overloaded division operator.
        """
        assert(self.__class__.__name__ == other.__class__.__name__)
        _x = numpy.union1d(self.x, other.x)
        _x = _x[other(_x) != 0]
        _y = self(_x)/other(_x)
        return self.__class__(_x, _y)

    def __add__(self, other):
        """Overloaded sum operator.
        """
        assert(self.__class__.__name__ == other.__class__.__name__)
        _x = numpy.union1d(self.x, other.x)
        _y = self(_x) + other(_x)
        return self.__class__(_x, _y)

    def __sub__(self, other):
        """Overloaded sum operator.
        """
        assert(self.__class__.__name__ == other.__class__.__name__)
        _x = numpy.union1d(self.x, other.x)
        _y = self(_x) - other(_x)
        return self.__class__(_x, _y)

    def __len__(self):
        """Return the lenght of the arrays used to construct the spline.
        """
        return len(self.x)

    def scale(self, scale, yname=None, yunits=None):
        """Scale the spline y values.

        Warning
        -------
        Need to set the names and units properly---the main issue is ths for
        derived classes the names of the corresponding members can be
        different.
        """
        _x = numpy.copy(self.x)
        _y = numpy.copy(self.y)*scale
        if yname is None:
            yname = self.yname
        if yunits is None:
            yunits = self.yunits
        fmt = dict(xname=self.xname, xunits=self.xunits, yname=yname,
                   yunits=yunits)
        return self.__class__(_x, _y, **fmt)

    @classmethod
    def label(cls, name, units=None):
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

        Note that the cdf is built using a linear interpolated spline, no matter
        what the class of the original spline is.
        """
        _x = self.x
        _y = numpy.array([self.integral(_x[0], _xp) for _xp in _x])/self.norm()
        return xInterpolatedUnivariateSplineLinear(_x, _y)

    def build_ppf(self):
        """Create the percent point function (or inverse of the cdf).

        Note that the cdf is built using a linear interpolated spline, no matter
        what the class of the original spline is.
        """
        _y = self.x
        _x = numpy.array([self.integral(_y[0], _xp) for _xp in _y])
        _x, _mask = numpy.unique(_x, return_index=True)
        _x/= self.norm()
        _y = _y[_mask]
        fmt = dict(xname='Normalized integral', yname=self.xname,
                   yunits=self.xunits)
        return xInterpolatedUnivariateSplineLinear(_x, _y, **fmt)

    def plot(self, num_points=1000, overlay=False, logx=False, logy=False,
             scale=1., offset=0., show=True, **kwargs):
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
        if not logx:
            _x = numpy.linspace(self.xmin(), self.xmax(), num_points)
        else:
            _x = numpy.logspace(numpy.log10(self.xmin()),
                                numpy.log10(self.xmax()), num_points)
        _y = scale*self(_x) + offset
        if overlay:
            plt.plot(_x, _y, '-', self.x, self.y, 'o', **kwargs)
        else:
            plt.plot(_x, _y, '-', **kwargs)
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

    This particular class is offering the facility to optimize the input
    arrays via the `optimize` argument.

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

    optimize : bool
        If `True`, the input arrays are optimized via the\
        `optimize_grid_linear()` function.

    tolerance : float
        The tolerance for the input array optimization. (If `optimize` is\
        `False`, this has no effect.)

    Example
    -------
    >>> from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
    >>>
    >>> x = numpy.linspace(0, 2*numpy.pi, 100)
    >>> y = numpy.sin(x)
    >>> s = xInterpolatedUnivariateSplineLinear(x, y, xname='x', xunits='au')
    >>> s.plot()
    """

    def __init__(self, x, y, xname=None, xunits=None, yname=None, yunits=None,
                 optimize=False, tolerance=1e-4):
        """ Constructor.
        """
        if optimize:
            oldx, oldy = x, y
            x, y = optimize_grid_linear(x, y, tolerance)
        xInterpolatedUnivariateSpline.__init__(self, x, y, None, [None, None],
                                               1, xname, xunits, yname, yunits)
        if optimize:
            dist = self.dist(oldx, oldy)
            logger.info('Relative (max/ave) dist. to original array: %e/%e' %\
                        (dist.max(), dist.sum()/len(dist)))


class xInterpolatedUnivariateLogSpline(xUnivariateSplineBase, UnivariateSpline):

    """
    """

    def __init__(self, x, y, w=None, bbox=[None, None], k=3,
                 xname=None, xunits=None, yname=None, yunits=None):
        """Constructor.
        """
        xUnivariateSplineBase.__init__(self, x, y, xname, xunits, yname, yunits)
        _x = numpy.log10(x)
        _y = numpy.log10(y)
        UnivariateSpline.__init__(self, _x, _y, w, bbox, k, s=None)
        
    def __call__(self, x):
        """Overloaded call method.
        """
        return numpy.power(10., UnivariateSpline.__call__(self, numpy.log10(x)))

    def integral(self, x1, x2):
        """Overloaded integral method.
        """
        from scipy.interpolate import dfitpack
        tck = self._eval_args
        return dfitpack.splint(*(tck + (x1, x2)))


class xInterpolatedUnivariateLogSplineLinear(xInterpolatedUnivariateLogSpline):

    """
    """

    def __init__(self, x, y, xname=None, xunits=None, yname=None, yunits=None):
        """Constructor.
        """
        xInterpolatedUnivariateLogSpline.__init__(self, x, y, None,
                                                  [None, None], 1, xname,
                                                  xunits, yname, yunits)


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

    def scale(self, scale_factor):
        """
        """
        return self.__class__(self.x.copy(), self.y.copy(), self.z*scale_factor)


class xInterpolatedBivariateSpline(xBivariateSplineBase, RectBivariateSpline):

    """Bivariate interpolated spline on a rectangular grid.
    """

    def __init__(self, x, y, z, kx=1, ky=1, xname=None, xunits=None,
                 yname=None, yunits=None, zname=None, zunits=None):
        """Constructor.
        """
        if hasattr(z, '__call__'):
            _x, _y = numpy.meshgrid(y, x)
            z = z(_x, _y)
        xBivariateSplineBase.__init__(self, x, y, z, xname, xunits, yname,
                                      yunits, zname, zunits)
        RectBivariateSpline.__init__(self, x, y, z,
                                     bbox=[None, None, None, None],
                                     kx=kx, ky=ky, s=0)

    def __call__(self, x, y, dx=0, dy=0, grid=False):
        """Overloaded __call__method.

        Here we basically override the default value of the `grid` parameter
        from `True` to `False`, since we're typically interested in evaluating
        the splined at given physical coordinates, rather than grid points.
        """
        return RectBivariateSpline.__call__(self, x, y, None, dx, dy, grid)

    def vslice(self, x):
        """Return a vertical slice at a given x of the bivariate spline.

        Args
        ----
        x : float
            The x value at which the vertical slice should be calculated.
        """
        _x = self.y
        _y = self(x, _x)
        fmt = dict(xname=self.yname, xunits=self.yunits, yname=self.zname,
                   yunits=self.zunits)
        return xInterpolatedUnivariateSplineLinear(_x, _y, **fmt)

    def hslice(self, y):
        """Return an horizontal slice at a given y of the bivariate spline.

        Args
        ----
        y : float
            The y value at which the horizontal slice should be calculated.
        """
        _x = self.x
        _y = self(_x, y)
        fmt = dict(xname=self.xname, xunits=self.xunits, yname=self.zname,
                   yunits=self.zunits)
        return xInterpolatedUnivariateSplineLinear(_x, _y, **fmt)

    def build_vppf(self):
        """Create the vertical percent point function (or inverse of cdf).

        Warning
        -------
        This really, really need to be fixed. Instead of grabbing a vertical
        slice at xmean, we should pass an argument to the function so that
        the subclasses can implement whatever is right for them.
        """
        _xmean = 0.5*(self.xmin() + self.xmax())
        _refppf = self.vslice(_xmean).build_ppf()
        _x = self.x.copy()
        _y = _refppf.x
        _z = numpy.zeros(shape = (_x.size, _y.size))
        for i, _xp in enumerate(_x):
            _ppf = self.vslice(_xp).build_ppf()
            for j, _yp in enumerate(_y):
                _z[i, j] = _ppf(_yp)
        fmt = dict(yname='Normalized integral', xname=self.xname,
                   xunits=self.xunits, zname=self.yname, zunits=self.yunits)
        return xInterpolatedBivariateSplineLinear(_x, _y, _z, **fmt)

    def plot(self, num_pointsx=100, num_pointsy=100, num_contours=75,
             logz=False, show=True):
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
        if logz:
            from matplotlib.colors import LogNorm
            _levels = numpy.logspace(-4, numpy.log10(_z.max()), num_contours)
            contour = plt.contourf(_x, _y, _z, levels=_levels, norm = LogNorm())
        else:
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

            
class xInterpolatedBivariateSplineLinear(xInterpolatedBivariateSpline):

    """Bivariate linear interpolated spline on a rectangular grid.
    """

    def __init__(self, x, y, z, xname=None, xunits=None,
                 yname=None, yunits=None, zname=None, zunits=None):
        xInterpolatedBivariateSpline.__init__(self, x, y, z, 1, 1, xname,
                                              xunits, yname, yunits, zname,
                                              zunits)


def main():
    x = numpy.linspace(0, 2*numpy.pi, 100)
    y = numpy.sin(x)
    s = xInterpolatedUnivariateSplineLinear(x, y, 'x', 'au', 'y')
    s.plot()


if __name__ == '__main__':
    main()
