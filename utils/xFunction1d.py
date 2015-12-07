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
import scipy.interpolate
import scipy.integrate

from ximpol.__logging__ import logger



class xFunction1d(scipy.interpolate.interp1d):

    """ Light-weight wrapper over the scipy scipy.interpolate.interp1d class.
    """

    def __init__(self, x, y, kind = 'quadratic'):
        """ Constructor.

        x and y are arrays of the tabulated function points and the kind
        argument specifies the type of interpolation as a string, the options
        being:
        - 'linear'   : linear interpolation
        - 'nearest'  : nearest value
        - 'zero'     : ?
        - 'slinear'  : first-order spline
        - 'quadratic': second-order spline
        - 'cubic'    : third-order spline

        Note that we explicitely set the assume_sorted argument of the
        scipy.interpolate.interp1d to false.
        """
        scipy.interpolate.interp1d.__init__(self, x, y, kind,
                                            assume_sorted = False)
        self.__Normalization = None

    def __mul__(self, other):
        """ Multiply two functions.
        """
        pass

    def __len__(self):
        """ Return the lenght of the underlying arrays.
        """
        return len(self.x)

    def optimize(self, rtol = 0.01, atol = 0):
        """ Optimize the sampling grid for a function.
        
        TODO: need to work on this one.
        """
        logger.info('Optimizing grid (rtol = %e, atol = %e)...' % (rtol, atol))
        newx = [self.x[0]]
        newy = [self.y[0]]
        for _x, _y in zip(self.x, self.y)[1:-1]:
            delta = abs(_y - newy[-1])
            ave = 0.5*(_y + newy[-1])
            if delta/ave > rtol or delta > atol:
                newx.append(_x)
                newy.append(_y)
        newx.append(self.x[-1])
        newy.append(self.y[-1])
        logger.info('Done, %d samples reduced to %d.' % (len(self), len(newx)))
        scipy.interpolate.interp1d.__init__(self, newx, newy,
                                            assume_sorted = False)

    def xmin(self):
        """ Return the minimun of the function support.

        (Since the x array is sorted, we just return its first element).
        """
        return self.x[0]

    def xmax(self):
        """ Return the maximum of the function support.
        
        (Since the x array is sorted, we just return its last element).
        """
        return self.x[-1]

    def evaluate(self, x):
        """ Convenience method to evaluate the function at a x.
        """
        return self(x)

    def integral(self, x1, x2):
        """ Evaluate the definite integral of the function over a given
        interval.
        """
        return scipy.integrate.quad(self, x1, x2)[0]

    def norm(self):
        """ Evaluate the normalization of the function, i.e., the integral
        over the entire support.

        Note that we cache the value of the integral the first time the
        method is called, in order to avoid making the calculation over and
        over again.
        """
        if self.__Normalization is None:
            self.__Normalization = self.integral(self.xmin(), self.xmax())
        return self.__Normalization

    def draw(self, npoints = 200):
        """ Draw the function.
        """
        import matplotlib.pyplot as plt
        _x = numpy.linspace(self.xmin(), self.xmax(), npoints)
        _y = self.evaluate(_x)
        plt.plot(_x, _y, '-', self.x, self.y, 'o')
        plt.show()

        

def test():
    """ Test code.
    """
    x = numpy.linspace(0, 2*numpy.pi, 20)
    y = numpy.sin(x)
    f1 = xFunction1d(x, y, 'linear')
    f2 = xFunction1d(x, y, 'quadratic')
    _x = numpy.array([0, numpy.radians(30), numpy.radians(45),
                      numpy.radians(60), 0.5*numpy.pi])
    _y = numpy.array([0, 0.5, 0.5*numpy.sqrt(2), 0.5*numpy.sqrt(3), 1])
    print(_y)
    print(f1(_x) - _y)
    print(f2(_x) - _y)
    print(f2.integral(0, numpy.pi))
    print(f2.norm())
    f2.draw()
    


if __name__ == '__main__':
    test()
