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


"""Unit test for the energy dispersion facility.
"""

import os
import numpy
import unittest
import sys

from ximpol import XIMPOL_IRF
from ximpol.irf.rmf import xEnergyDispersion
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.matplotlib_ import overlay_tag, save_current_figure
from ximpol.detector.xipe import _full_path
from ximpol.irf import load_rmf


IRF_NAME = 'xipe_baseline'
GPD_ERES_FILE_PATH = _full_path('eres_fwhm_hedme8020_1atm_1cm.asc')



class TestEnergyDispersion(unittest.TestCase):

    """Unit test for xModulationFactor.
    """

    @classmethod
    def setUpClass(self):
        """Setup---here we essentially create the modulation factor.
        """
        self.measx, self.measy = numpy.loadtxt(GPD_ERES_FILE_PATH, unpack=True)
        self.edisp = load_rmf(IRF_NAME)
        self.interactive = sys.flags.interactive

    def test_rvs(self, num_events=100000):
        """
        """
        mc_energy = 10
        _ppf = self.edisp.matrix.vppf.vslice(mc_energy)
        _e = numpy.zeros(num_events)
        _e.fill(mc_energy)
        _slice = self.edisp.matrix.slice(mc_energy)
        _ch = self.edisp.matrix.rvs(_e)
        n, bins, patches = plt.hist(_ch, bins=numpy.linspace(0, 255, 256),
                                    histtype='step')
        bin_width = (bins[1] - bins[0])
        scale = n.sum()*bin_width/_slice.norm()
        plt.plot(_slice.x, scale*_slice.y)
        #plt.show()


if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)
