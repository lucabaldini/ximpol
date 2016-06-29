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

from ximpol.irf.mrf import xAzimuthalResponseGenerator
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.irf import load_arf, load_mrf, DEFAULT_IRF_NAME

"""
"""

class TestPolarization(unittest.TestCase):

    """
    """

    @classmethod
    def setUpClass(cls):
        """Setup.
        """
        cls.aeff = load_arf(DEFAULT_IRF_NAME)
        cls.modf = load_mrf(DEFAULT_IRF_NAME)
        cls.generator = xAzimuthalResponseGenerator()
        cls.emin = 1.
        cls.emax = 10.

    def polarization_degree(self, energy):
        """
        """
        return (1. - (energy - self.emin)/(self.emax - self.emin))

    def test_simplest(self, size=100000, phase=0.):
        """
        """
        modf = 1.
        energy = numpy.random.uniform(self.emin, self.emax, size)
        visibility = modf*self.polarization_degree(energy)
        phi = self.generator.rvs_phi(visibility, phase)
        binning = numpy.linspace(0, 2*numpy.pi, 100)
        hist = plt.hist(phi, bins=binning)
        fit_results = self.generator.fit_histogram(hist)
        fit_results.set_polarization(modf)
        fit_degree = fit_results.polarization_degree
        mean_energy = numpy.mean(energy)
        exp_degree = self.polarization_degree(mean_energy)
        delta = abs(fit_degree - exp_degree)/exp_degree
        msg = 'delta = %.3f' % delta
        self.assertTrue(delta < 0.05, msg)

    def test_uniform(self, size=100000, phase=0.):
        """
        """
        energy = numpy.random.uniform(self.emin, self.emax, size)
        visibility = self.modf(energy)*self.polarization_degree(energy)
        phi = self.generator.rvs_phi(visibility, phase)
        binning = numpy.linspace(0, 2*numpy.pi, 100)
        hist = plt.hist(phi, bins=binning)
        fit_results = self.generator.fit_histogram(hist)
        mean_energy = numpy.mean(energy)
        mu_effective = self.modf.weighted_average(energy)
        fit_results.set_polarization(mu_effective)
        fit_degree = fit_results.polarization_degree
        exp_degree = (self.polarization_degree(energy)*self.modf(energy))\
                     .sum()/self.modf(energy).sum()
        delta = abs(fit_degree - exp_degree)/exp_degree
        msg = 'delta = %.3f' % delta
        self.assertTrue(delta < 0.05, msg)



if __name__ == '__main__':
    unittest.main()
