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


"""Unit test for the count spectrum facility.
"""

import os
import numpy
import unittest

from ximpol import XIMPOL_IRF
from ximpol.test import save_current_figure
from ximpol.irf.arf import xEffectiveArea
from ximpol.srcmodel.spectrum import xCountSpectrum
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.utils.matplotlib_ import pyplot as plt


class TestCountSpectrum(unittest.TestCase):

    """Unit test for xCountSpectrum.
    """

    @classmethod
    def setUpClass(self):
        """Setup.
        """
        file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.arf')
        self.aeff = xEffectiveArea(file_path)

    def test_power_law_stationary(self, interactive=False):
        """Test a time-independent power law.

        This creates a count spectrum witn no time dependence, i.e., with
        only two (identical) interpolating time points on the auxiliary
        (time) axis. The power-law parameters (C and gamma) are constant and the
        underlying xUnivariateAuxGenerator looks like.

        .. image:: ../../ximpol/test/figures/test_power_law_stationary_2d.png

        Then a vertical slice (i.e., an interpolated linear spline) is taken
        in the middle of the auxiliary axis and the y-values of the spline
        are compared with the direct product of the effective area and the
        count spectrum. If everything goes well, they should be on top
        of each other. The figure below is also showing the orignal power-law
        spectrum multiplied by the peak effective area.

        .. image:: ../../ximpol/test/figures/test_power_law_stationary_slice.png

        """
        tmin = 0.
        tmax = 100.
        tref = 0.5*(tmin + tmax)
        C = 1.
        Gamma = 2.

        def powerlaw(E, t):
            """Function defining a time-dependent energy spectrum.
            """
            return C*numpy.power(E, -Gamma)

        _t = numpy.linspace(tmin, tmax, 2)
        count_spectrum = xCountSpectrum(powerlaw, self.aeff, _t)
        count_spectrum.plot(show=interactive)
        save_current_figure('test_power_law_stationary_2d.png')

        ref_slice = count_spectrum.slice(tref)
        _x = self.aeff.x
        _y = self.aeff.y.max()*C*numpy.power(_x, -Gamma)
        plt.plot(_x, _y, '-')
        _y = C*numpy.power(_x, -Gamma)*self.aeff.y
        plt.plot(_x, _y, 'o')
        ref_slice.plot(logx=True, logy=True, show=interactive)
        save_current_figure('test_power_law_stationary_slice.png')

        delta = abs((_y - ref_slice(_x))/_y).max()
        self.assertTrue(delta < 1e-3, 'max deviation %.9f' % delta)

    def test_power_law_variable(self, interactive=False):
        """Test a time-dependent power law.
        """
        tmin = 0.
        tmax = 100.
        tref = 0.5*(tmin + tmax)

        def C(t):
            """Time-dependent C---equals to 1 @ tmin and 2 @ tmax.
            """
            return 1. + (t - tmin)/(tmax - tmin)

        def Gamma(t):
            """Time-dependent C---equals to 2 @ tmin and 1 @ tmax.
            """
            return 2. - (t - tmin)/(tmax - tmin)

        def powerlaw(E, t):
            """Function defining a time-dependent energy spectrum.
            """
            return C(t)*numpy.power(E, -Gamma(t))


        _t = numpy.linspace(tmin, tmax, 100)
        count_spectrum = xCountSpectrum(powerlaw, self.aeff, _t)
        count_spectrum.plot(show=interactive)
        save_current_figure('test_power_law_variable_2d.png')

        ref_slice = count_spectrum.slice(tref)
        _x = self.aeff.x
        _y = self.aeff.y.max()*C(tref)*numpy.power(_x, -Gamma(tref))
        plt.plot(_x, _y, '-')
        _y = C(tref)*numpy.power(_x, -Gamma(tref))*self.aeff.y
        plt.plot(_x, _y, 'o')
        ref_slice.plot(logx=True, logy=True, show=interactive)
        save_current_figure('test_power_law_variable_slice.png')

        delta = abs((_y - ref_slice(_x))/_y).max()
        self.assertTrue(delta < 1e-3, 'max deviation %.9f' % delta)




if __name__ == '__main__':
    unittest.main()
