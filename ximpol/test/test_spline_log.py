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


import sys
import unittest

from ximpol.core.spline import *
from ximpol.utils.logging_ import suppress_logging
from ximpol.utils.matplotlib_ import pyplot as plt


class TestSplineLog(unittest.TestCase):

    """Unit test for xInterpolatedUnivariateSplineLinear.
    """

    @classmethod
    def setUpClass(cls):
        """Setup.
        """
        cls.interactive = sys.flags.interactive

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
        slin = xInterpolatedUnivariateSplineLinear(_x, _y)
        slog = xInterpolatedUnivariateLogSplineLinear(_x, _y)
        target_norm = self.power_law_integral(norm, index, emin, emax)
        lin_norm = slin.norm()
        log_norm = slog.norm()
        delta = abs(target_norm - log_norm)/target_norm
        msg = 'delta = %.3e' % delta
        self.assertTrue(delta < 0.01, msg)
        if self.interactive:
            plt.figure()
            slin.plot(logx=True, logy=True, overlay=True, show=False)
            slog.plot(logx=True, logy=True, overlay=True, show=False)



if __name__ == '__main__':
    unittest.main()
