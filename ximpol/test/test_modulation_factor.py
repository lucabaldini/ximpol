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


"""Unit test for the modulation factor facility.
"""

import os
import numpy
import unittest
import sys

from ximpol import XIMPOL_IRF
from ximpol.irf.mrf import xModulationFactor
from ximpol.irf.mrf import xAzimuthalResponseGenerator
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.matplotlib_ import overlay_tag, save_current_figure
from ximpol.detector.xipe import _full_path
from ximpol.irf import load_mrf


IRF_NAME = 'xipe_baseline'
GPD_MODF_FILE_PATH = _full_path('modfact_hedme8020_1atm_1cm_mng.asc')


"""We explictely set the random seed to have reproducible results.
"""
numpy.random.seed(0)


class TestModulationFactor(unittest.TestCase):

    """Unit test for xModulationFactor.
    """

    @classmethod
    def setUpClass(self):
        """Setup---here we essentially create the modulation factor.
        """
        self.measx, self.measy = numpy.loadtxt(GPD_MODF_FILE_PATH, unpack=True)
        self.modf = load_mrf(IRF_NAME)
        self.interactive = sys.flags.interactive

    def test_constant(self, num_events=1000000, polarization_degree=1.,
                      polarization_angle=numpy.radians(20)):
        """Test the modulation factor as a random number generator when
        both the polarization angle and degrees are energy- and
        time-independent.
        """
        poldegree = numpy.zeros(num_events)
        poldegree.fill(polarization_degree)
        polangle = numpy.zeros(num_events)
        polangle.fill(polarization_angle)
        self.modf.generator.plot(show=False)
        save_current_figure('test_modulation_constant_generator.png',
                            show=self.interactive)
        emin = self.modf.xmin()
        emax = self.modf.xmax()
        energy = numpy.random.uniform(emin, emax, num_events)
        phi = self.modf.rvs_phi(energy, poldegree, polangle)
        ebinning = numpy.linspace(emin, emax, 10)
        phi_binning = numpy.linspace(0, 2*numpy.pi, 100)
        fit_results = []
        for i, (_emin, _emax) in enumerate(zip(ebinning[:-1], ebinning[1:])):
            _emean = 0.5*(_emin + _emax)
            _mask = (energy > _emin)*(energy < _emax)
            _phi = phi[_mask]
            _hist = plt.hist(_phi, bins=phi_binning, histtype='step')
            _fr = xAzimuthalResponseGenerator.fit_histogram(_hist)
            _fr.emean = _emean
            fit_results.append(_fr)
            _fr.plot(label='Energy: %.2f--%.2f keV' % (_emin, _emax))
            plt.axis([0., 2*numpy.pi, 0., 1.2*_hist[0].max()])
            overlay_tag()
            save_current_figure('test_modulation_constant_fit_slice%d.png' % i,
                                show=self.interactive)
        _x = [_fr.emean for _fr in fit_results]
        _y = [_fr.phase for _fr in fit_results]
        _dy = [_fr.phase_error for _fr in fit_results]
        plt.errorbar(_x, _y, yerr=_dy, fmt='o')
        plt.plot(_x, numpy.array([polarization_angle]*len(_x)))
        plt.xlabel('Energy [keV]')
        plt.ylabel('Modulation angle [$^\circ$]')
        save_current_figure('test_modulation_constant_angle.png',
                            show=self.interactive)
        _y = [_fr.visibility for _fr in fit_results]
        _dy = [_fr.visibility_error for _fr in fit_results]
        plt.errorbar(_x, _y, yerr=_dy, fmt='o')
        plt.axis([emin, emax, 0, 1])
        self.modf.plot(show=False)
        plt.xlabel('Energy [keV]')
        plt.ylabel('Modulation visibility')
        save_current_figure('test_modulation_constant_visibility.png',
                            show=self.interactive)


if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)
