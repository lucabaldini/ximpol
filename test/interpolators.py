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


import unittest

from ximpol.core.xInterpolatedUnivariateSpline import *



class testInterpolatedUnivariateSplineLinear(unittest.TestCase):

    """ Unit test for xInterpolatedUnivariateSplineLinear.
    """

    def setUp(self):
        """Setup.

        Create a few objects to be used for testing.
        """
        self.num_points = 100
        self.x1 = numpy.linspace(0, 2*numpy.pi, self.num_points)
        self.y1 = numpy.sin(self.x1)
        self.x2 = numpy.linspace(0, numpy.pi, self.num_points)
        self.y2 = numpy.sin(self.x2)
        self.s1 = xInterpolatedUnivariateSplineLinear(self.x1, self.y1)
        self.s2 = xInterpolatedUnivariateSplineLinear(self.x2, self.y2)

    def test_basic(self):
        """Test the basic object instantiation.
        """
        # Check we get the number of points right.
        self.assertEqual(len(self.s1), self.num_points)

    def test_evaluation(self):
        """Test the object evaluation.
        """
        # This is a linear interpolator, so the interpolated values must
        # be identical, within rounding errors, to the original grid of
        # values.
        _delta = abs(self.s1(self.x1) - self.y1)
        self.assertTrue(_delta.all() < 1e-9)

        # s1 and s2 are built with different sets of points, but with the same
        # underlying function, so they should be fairly close at any given
        # point.
        _x = numpy.linspace(0, numpy.pi, 10)
        _delta = abs(self.s1(_x) - self.s2(_x))
        self.assertTrue(_delta.all() < 1e-6)

    def test_multiplication(self):
        """Test the interpolator multiplication.
        """
        # Evaluate s1*s2 in x2 should give the same answer than multiplying
        # s1(x2)*y2.
        _m = self.s1*self.s2
        _delta = abs(_m(self.x2) - self.s1(self.x2)*self.y2)
        self.assertTrue(_delta.all() < 1e-9)

        # And the result of the multiplication should be an instance of
        # the original operand class.
        self.assertTrue(isinstance(_m, xInterpolatedUnivariateSplineLinear))

    def test_sum(self):
        """Test the interpolator sum.
        """
        # Evaluate s1 + s2 in x2 should give the same answer than adding
        # s1(x2) + y2.
        _s = self.s1 + self.s2
        _delta = abs(_s(self.x2) - (self.s1(self.x2) + self.y2))
        self.assertTrue(_delta.all() < 1e-9)

        # And the result of the multiplication should be an instance of
        # the original operand class.
        self.assertTrue(isinstance(_s, xInterpolatedUnivariateSplineLinear))

    def test_extrapolation(self):
        """Test interpolator extrapolation.
        """
        pass



if __name__ == '__main__':
    unittest.main()
