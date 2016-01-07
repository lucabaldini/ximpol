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

from ximpol.core.spline import *
from ximpol.utils.logging_ import suppress_logging
suppress_logging()


class TestSplineLinear(unittest.TestCase):

    """Unit test for xInterpolatedUnivariateSplineLinear.
    """

    @classmethod
    def setUpClass(self):
        """Setup.

        Create a few objects to be used for testing.
        """
        self.num_points = 100
        self.x1 = numpy.linspace(0, 2*numpy.pi, self.num_points)
        self.y1 = numpy.sin(self.x1)
        self.x2 = numpy.linspace(0, numpy.pi, self.num_points)
        self.y2 = numpy.sin(self.x2)
        self.x3 = numpy.linspace(0, 10, 100)
        self.y3 = 3*self.x3
        self.s1 = xInterpolatedUnivariateSplineLinear(self.x1, self.y1)
        self.s2 = xInterpolatedUnivariateSplineLinear(self.x2, self.y2)
        self.s3 = xInterpolatedUnivariateSplineLinear(self.x3, self.y3)

    def test_len(self):
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
        self.assertTrue(_delta.max() < 1e-9, 'max. diff. %.9f' % _delta.max())

        # s1 and s2 are built with different sets of points, but with the same
        # underlying function, so they should be fairly close at any given
        # point.
        _x = numpy.linspace(0, numpy.pi, 10)
        _delta = abs(self.s1(_x) - self.s2(_x))
        self.assertTrue(_delta.max() < 1e-3, 'max. diff. %.9f' % _delta.max())

    def test_multiplication(self):
        """Test the interpolator multiplication.
        """
        # Evaluate s1*s2 in x2 should give the same answer than multiplying
        # s1(x2)*y2.
        _m = self.s1*self.s2
        _delta = abs(_m(self.x2) - self.s1(self.x2)*self.y2)
        self.assertTrue(_delta.max() < 1e-9, 'max. diff. %.9f' % _delta.max())

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
        self.assertTrue(_delta.max() < 1e-9, 'max. diff. %.9f' % _delta.max())

        # And the result of the multiplication should be an instance of
        # the original operand class.
        self.assertTrue(isinstance(_s, xInterpolatedUnivariateSplineLinear))

    def test_extrapolation(self):
        """Test interpolator extrapolation.
        """
        # Calculate one extrapolated value by hand and compare it to the
        # value from the interpolator.
        _xa = self.x1[-2]
        _xb = self.x1[-1]
        _ya = self.y1[-2]
        _yb = self.y1[-1]
        _x = _xb + 0.2
        _y = _ya + (_yb - _ya)/(_xb - _xa)*(_x - _xa)
        _delta = abs(self.s1(_x) - _y)
        self.assertTrue(_delta < 1e-9, 'max. diff. %.9f' % _delta)

    def test_norm(self):
        """Test the normalization calculation.
        """
        _delta = abs(self.s3.norm() - 100.*3./2)
        self.assertTrue(_delta < 1e6, 'norm. diff. %.9f' % _delta)

    def test_cdf(self):
        """ The cdf must be 0 at xmin and 1 at xmax.
        """
        cdf = self.s3.build_cdf()
        _delta = abs(cdf(self.s3.xmin()))
        self.assertTrue(_delta < 1e-3, 'ppf(xmin) %.9f' % _delta)
        _delta = abs(cdf(self.s3.xmax()) - 1.)
        self.assertTrue(_delta < 1e-3, 'ppf(xmax) - 1 %.9f' % _delta)

    def test_ppf(self):
        """ The ppf must be defined between 0 and 1 (where is equal to the
        xmin and xmax values of the original spline).
        """
        ppf = self.s3.build_ppf()
        _delta = abs(ppf.xmin())
        self.assertTrue(_delta < 1e-3, 'ppf xmin %.9f' % _delta)
        _delta = abs(ppf.xmax() - 1.)
        self.assertTrue(_delta < 1e-3, 'ppf (xmax - 1) %.9f' % _delta)
        _delta = abs(ppf(0) - self.s3.xmin())
        self.assertTrue(_delta < 1e-3, 'ppf(0) - xmin %.9f' % _delta)
        _delta = abs(ppf(1) - self.s3.xmax())
        self.assertTrue(_delta < 1e-3, 'ppf(1) - xmax %.9f' % _delta)


if __name__ == '__main__':
    unittest.main()
