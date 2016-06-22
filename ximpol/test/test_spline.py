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
suppress_logging()


class TestSplineLinear(unittest.TestCase):

    """Unit test for xInterpolatedUnivariateSplineLinear.
    """

    @classmethod
    def setUpClass(cls):
        """Setup.

        Create a few objects to be used for testing.
        """
        cls.num_points = 100
        cls.x1 = numpy.linspace(0, 2*numpy.pi, cls.num_points)
        cls.y1 = numpy.sin(cls.x1)
        cls.x2 = numpy.linspace(0, numpy.pi, cls.num_points)
        cls.y2 = numpy.sin(cls.x2)
        cls.x3 = numpy.linspace(0, 10, 100)
        cls.y3 = 3*cls.x3
        cls.s1 = xInterpolatedUnivariateSplineLinear(cls.x1, cls.y1)
        cls.s2 = xInterpolatedUnivariateSplineLinear(cls.x2, cls.y2)
        cls.s3 = xInterpolatedUnivariateSplineLinear(cls.x3, cls.y3)

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

    def test_cdf_erf(self):
        """Test the cdf for a gaussian function.
        """
        _x = numpy.linspace(-5, 5, 100)
        _y = 1./numpy.sqrt(2.*numpy.pi)*numpy.exp(-0.5*_x**2)
        pdf = xInterpolatedUnivariateSplineLinear(_x, _y)
        cdf = pdf.build_cdf()
        delta = abs(cdf(_x) - 0.5*(1. + scipy.special.erf(_x/numpy.sqrt(2.))))
        max_delta = delta.max()
        err_msg = 'maximum absolute delta %.4e' % max_delta
        self.assertTrue(max_delta < 5e-4, err_msg)

    def test_sort(self):
        """Test the automatic sorting functionality.
        """
        _x = numpy.random.sample(100)
        _y = _x**2
        s = xInterpolatedUnivariateSplineLinear(_x, _y)
        _x.sort()
        self.assertTrue((s.x == _x).all())
        self.assertTrue((s.y == _x**2).all())

    def test_non_unique(self):
        """The spline constructor must fail when non-unique values are passed.
        """
        _x = numpy.array([1, 1, 2, 3, 4])
        _y = _x**2
        with self.assertRaises(AssertionError):
            s = xInterpolatedUnivariateSplineLinear(_x, _y)
        
        

if __name__ == '__main__':
    unittest.main()
