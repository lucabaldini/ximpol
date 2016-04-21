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


"""Unit test for the irf.rmf module.
"""

import unittest
import os
import numpy

from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.detector.xipe import _full_path
from ximpol.irf import load_rmf


IRF_NAME = 'xipe_baseline'
GPD_ERES_FILE_PATH = _full_path('eres_fwhm_hedme8020_1atm_1cm.asc')


class TestXipeRmf(unittest.TestCase):

    """Unit test for the XIPE energy dispersion.
    """

    @classmethod
    def setUpClass(cls):
        """Setup.

        Create a few objects to be used for testing.
        """
        cls.edisp = load_rmf(IRF_NAME)

    def test_xipe_rmf_matrix_norm(self):
        """Test the XIPE energy dispersion normalization.

        Take a number of vertical slices of the energy dispersion matrix
        and make sure they are normalized to unity.
        """
        emin = self.edisp.matrix.xmin()
        emax = self.edisp.matrix.xmax()
        de = emax - emin
        for energy in numpy.linspace(emin + 0.2*de, emax - 0.2*de, 10):
            _delta = abs(self.edisp.matrix.vslice(energy).norm() - 1)
            self.assertTrue(_delta < 1e-3, 'diff. %.9f' % _delta)

    def test_xipe_rmf_matrix_sigma(self):
        """Test the XIPE energy resolution.

        Load the text file with the XIPE FWHM energy resolution,
        recalculate the FWHM from the actual energy dispersion, as read from
        the corresponding .rmf FITS file, and make sure the two are
        sufficiently close to each other.
        """
        _x, _y = numpy.loadtxt(GPD_ERES_FILE_PATH, unpack=True)
        for energy, fwhm in zip(_x, _y):
            _slice = self.edisp.matrix.vslice(energy)
            _ppf = _slice.build_ppf()
            _sigma = 0.5*(self.edisp.ebounds(_ppf(0.8413)) -\
                          self.edisp.ebounds(_ppf(0.1586)))
            _delta = abs((fwhm - 2.358*_sigma)/fwhm)
            self.assertTrue(_delta < 2e-2, 'diff. %.9f' % _delta)


if __name__ == '__main__':
    unittest.main()
