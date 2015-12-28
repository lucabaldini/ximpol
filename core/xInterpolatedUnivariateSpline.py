#!/usr/bin/env python
# *********************************************************************
# * Copyright (C) 2015 Luca Baldini (luca.baldini@pi.infn.it)         *
# *                                                                   *
# * For the license terms see the file LICENSE, distributed           *
# * along with this software.                                         *
# *********************************************************************
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
from scipy.interpolate import InterpolatedUnivariateSpline

from ximpol.__logging__ import logger



class xInterpolatedUnivariateSpline(InterpolatedUnivariateSpline):

    """Light-weight wrapper over the standard InterpolatedUnivariateSpline
    class provided by scipy.interpolate, see
    http://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.InterpolatedUnivariateSpline.html

    The basic additional features we are implementing here are:
    1\ we do keep track of the original arrays passed to the interpolator;
    2\ optional xmin and xmax parameters can be passed to the constructor
       to limit the range over which the interpolator is defined;
    3\ multiplication is supported.

    Note that the interface to the base class has changed from numpy 0.14.
    An ext argument can be passed to the constructor starting with scipy
    0.15 to control the extrapolation behavior and a check_finite argument is
    available in 0.16. We currently do not use either one.
    """

    def __init__(self, x, y = None, w=None, bbox=[None, None], k=3,
                 xmin=-numpy.inf, xmax=numpy.inf):
        """Constructor.

        The arguments are the same of those taken by the native scipy class,
        except for xmin and xmax.

        If x is a string, it is interpreted as a path to a txt file name
        from which the data points are loaded.
        """
        if isinstance(x, str):
            file_path = x
            logger.info('Reading data values from %s...' % file_path)
            x, y = numpy.loadtxt(file_path, unpack = True)
        assert(len(x) == len(y))
        _mask = (x >= xmin)*(x <= xmax)
        self.__x, self.__y = x[_mask], y[_mask]
        InterpolatedUnivariateSpline.__init__(self, x, y, w, bbox, k)

    def __xmerge(self, other):
        """
        """
        _xmin = max(self.xmin(), other.xmin())
        _xmax = min(self.xmax(), other.xmax())
        assert(_xmax > _xmin)
        _x = numpy.concatenate((self.x(), other.x()))
        _x.sort()
        return _x

    def __mul__(self, other):
        """Multiplication operator.
        """
        assert(self.__class__.__name__ == other.__class__.__name__)
        _x = self.__xmerge(other)
        _y = self(_x)*other(_x)
        return self.__class__(_x, _y)

    def __add__(self, other):
        """Addition operator.
        """
        assert(self.__class__.__name__ == other.__class__.__name__)
        _x = self.__xmerge(other)
        _y = self(_x) + other(_x)
        return self.__class__(_x, _y)

    def __len__(self):
        """Return the lenght of the underlying arrays used to construct the
        interpolator.
        """
        return len(self.__x)

    def x(self):
        """Return the x array used to construct the interpolator.
        """
        return self.__x

    def y(self):
        """Return the y array used to construct the interpolator.
        """
        return self.__y

    def xmin(self):
        """Return the minimun of the x array used to construct the
        interpolator.
        """
        return self.__x[0]

    def xmax(self):
        """Return the maximum of the x array used to construct the
        interpolator..
        """
        return self.__x[-1]

    def plot(self, npts=1000):
        """Plot the function.
        """
        import matplotlib.pyplot as plt
        _x = numpy.linspace(self.xmin(), self.xmax(), npts)
        _y = self(_x)
        plt.plot(_x, _y, '-', self.__x, self.__y, 'o')
        plt.show()



class xInterpolatedUnivariateSplineLinear(xInterpolatedUnivariateSpline):

    """xInterpolatedUnivariateSpline subclass implementing the simplest
    poassible linear interpolator.
    """
    def __init__(self, x, y = None, xmin=-numpy.inf, xmax=numpy.inf):
        """ Constructor.
        """
        xInterpolatedUnivariateSpline.__init__(self, x, y, None, [None, None],
                                               1, xmin, xmax)



def test():
    """ Test code.
    """
    x = numpy.linspace(0, 2*numpy.pi, 20)
    y = numpy.sin(x)
    f1 = xInterpolatedUnivariateSplineLinear(x, y)
    f1.plot()



if __name__ == '__main__':
    test()
