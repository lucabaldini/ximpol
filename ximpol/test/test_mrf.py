#!/usr/bin/env python
#
# Copyright (C) 2016, the ximpol team.
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


"""Unit test for the irf.mrf module.
"""

import unittest
import os
import numpy

from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.detector.xipe import _full_path
from ximpol.irf import load_mrf


IRF_NAME = 'xipe_baseline'
GPD_MODF_FILE_PATH = _full_path('modfact_hedme8020_1atm_1cm_mng.asc')


class TestXipeMrf(unittest.TestCase):

    """Unit test for the Xipe modulation factor.
    """

    def test_xipe_mrf(self):
        """Test the XIPE modulation factor.

        This loads the modulation factor from the XIPE .mrf FITS file,
        then the corresponding values from the text file that the response
        functions are created from and makes sure that the two are close
        enough.
        """
        _x, _y = numpy.loadtxt(GPD_MODF_FILE_PATH, unpack=True)
        modf = load_mrf(IRF_NAME)
        _mask = _y > 0.
        _x = _x[_mask]
        _y = _y[_mask]
        _delta = abs((_y - modf(_x))/_y)
        self.assertTrue(_delta.max() < 5e-3, 'max. diff. %.9f' % _delta.max())


if __name__ == '__main__':
    unittest.main()
