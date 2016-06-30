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


"""Unit test for the core.spline module.
"""


import unittest
import scipy.special

from ximpol.core.spline import *
from ximpol.utils.logging_ import suppress_logging



class TestSplineLog(unittest.TestCase):

    """Unit test for xInterpolatedUnivariateSplineLinear.
    """

    @classmethod
    def setUpClass(cls):
        """Setup.
        """
        pass

    def power_law_integral(self, norm, index, xmin, xmax):
        """
        """
        return norm/(1. - index)*(xmax**(1. -index) - xmin**(1. - index))

    def test_power_law(self):
        """
        """
        norm = 1.
        index = 2.
        emin = 1.
        emax = 10.
        num_points = 5
        _x = numpy.logspace(numpy.log10(emin), numpy.log10(emax), num_points)
        _y = norm*_x**(-index)
        s = xUnivariateSplineLogLog(_x, _y)
        print self.power_law_integral(norm, index, emin, emax)
        print s.norm()
        s.plot(logx=True, logy=True, overlay=True)
        



if __name__ == '__main__':
    unittest.main()
