#!/urs/bin/env python
#
# Copyright (C) 2015--2016, the ximpol team.
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


import sys
import os
import unittest
import numpy

from ximpol.irf.mrf import xStokesAccumulator
from ximpol.irf.mrf import xAzimuthalResponseGenerator



class TestStokesAccumulator(unittest.TestCase):

    """Unit test for the Stokes accumulator.
    """

    def test_accumulator(self):
        """
        """
        accumulator = xStokesAccumulator()
        generator = xAzimuthalResponseGenerator()
        visibility = numpy.full(1000000, 0.43)
        phase = numpy.radians(45.)
        print visibility[0], phase
        phi = generator.rvs_phi(visibility, phase)
        accumulator.fill(phi)
        _vis, _vis_err = accumulator.visibility()
        _pha, _pha_err = accumulator.phase()
        r_vis = (_vis - visibility[0])/_vis_err
        r_pha = (_pha - phase)/_pha_err
        self.assertTrue(abs(r_vis) < 5., 'visibility residual %s' % r_vis)
        self.assertTrue(abs(r_pha) < 5., 'phase residual %s' % r_pha)
 
    
if __name__ == '__main__':
    unittest.main()
