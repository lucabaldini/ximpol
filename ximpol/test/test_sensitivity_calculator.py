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

from ximpol.srcmodel.spectrum import xCountSpectrum, power_law
from ximpol.srcmodel.spectrum import int_eflux2pl_norm
from ximpol.utils.units_ import keV2erg
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.irf import load_arf, load_mrf


"""The calculator ouputs have been obtained running by hand the code on the web
http://www.isdc.unige.ch/xipe/index.php/sensitivity-calculator
and are stored in the form of a list of tuples containing in order:
1 - column_density [1e22 cm^{-2}]
2 - power-law spectral index
3 - exposure time [ks]
4 - integral energy flux between 2 and 8 keV [1e-8 erg/cm^{2}/s]
5 - MDP in the 2--4, 4--6 and 6--8 keV energy bands.

Note that the numbers in the tuple are exactly what you would write in the
web form.
"""
SENSITIVITY_CALCULATOR_OUTPUT = [
    (0.1, 1., 10., 0.1, [0.04022, 0.06668, 0.14058]),
    (0.1, 2., 10., 0.1, [0.03293, 0.06927, 0.17443])
]


class TestSensitivityCalculator(unittest.TestCase):

    """Unit test for Fabio's sensitivity calculator at
    http://www.isdc.unige.ch/xipe/index.php/sensitivity-calculator
    """

    @classmethod
    def setUpClass(cls):
        """Setup.
        """
        cls.irf_name = 'xipe_baseline'
        cls.aeff = load_arf(cls.irf_name)
        cls.modf = load_mrf(cls.irf_name)
        cls.emin = 2.
        cls.emax = 8.
        cls.ebinning = numpy.linspace(cls.emin, cls.emax, 4)

    def mdp_table(self, column_density, index, exposure_time, flux):
        """Return the MDP table for a point source with a power-law
        spectral shape with a given set of parameters and for a given
        observation time.
        """
        norm = int_eflux2pl_norm(flux, self.emin, self.emax, index, erg=True)
        energy_spectrum = power_law(norm, index)
        tsamples = numpy.linspace(0., exposure_time, 2)
        count_spectrum = xCountSpectrum(energy_spectrum, self.aeff, tsamples)
        mdp_table = count_spectrum.build_mdp_table(self.ebinning, self.modf)
        return mdp_table

    def test_mdp(self):
        """Compare the MDP calculated by ximpol with that returned by the
        sensitivity calculator.
        """
        for column_density, index, exposure_time, flux, target_mdps in\
            SENSITIVITY_CALCULATOR_OUTPUT:
            # Convert the source parameters to physical units.
            column_density *= 1.e22
            exposure_time *= 1000.
            flux *= 1.e-8
            # Calculate the MDP table using the ximpol facilities.
            table = self.mdp_table(column_density, index, exposure_time, flux)
            ximpol_mdps = table.mdp_values()[:-1]
            target_mdps = numpy.array(target_mdps)
            ximpol_mdps = numpy.array(ximpol_mdps)
            delta = abs(target_mdps - ximpol_mdps)/target_mdps
            max_delta = delta.max()
            self.assertTrue(max_delta < 0.05, 'max. diff. %.4f' % max_delta)
        


if __name__ == '__main__':
    unittest.main()

