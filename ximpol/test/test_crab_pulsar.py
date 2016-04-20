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


"""Unit test for the Crab pulsar model.
"""

import os
import numpy
import unittest
import sys

from ximpol.core.rand import xUnivariateAuxGenerator
from ximpol.config.crab_pulsar import energy_spectrum, pl_normalization_spline,\
    pl_index_spline
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.matplotlib_ import overlay_tag, save_current_figure


"""We explictely set the random seed to have reproducible results.
"""
numpy.random.seed(0)


class TestCrabPulsarSpectrum(unittest.TestCase):

    """Unit test for xCountSpectrum.
    """

    @classmethod
    def setUpClass(cls):
        """Setup method.
        """
        pass

    def test_source_spectrum(self):
        """Test the intrinsic source spectrum.
        """
        energy = numpy.linspace(1., 10., 100)
        phase = numpy.linspace(0., 1., 100)
        fmt = dict(rvname='Energy', rvunits='keV', auxname='Phase')
        spec = xUnivariateAuxGenerator(phase, energy, energy_spectrum, **fmt)
        for _phi in numpy.linspace(0., 1., 33):
            _C = pl_normalization_spline(_phi)
            _Gamma = pl_index_spline(_phi)
            _pl = _C*numpy.power(energy, -_Gamma)
            _slice = spec.slice(_phi)(energy)
            delta = abs((_pl - _slice)/_pl).max()
            self.assertTrue(delta < 2.e-2, 'max deviation %.9f' % delta)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)
