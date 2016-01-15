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

from ximpol import XIMPOL_IRF
from ximpol.detector.xipe import GPD_MODF_FILE_PATH, IRF_NAME
from ximpol.irf.mrf import xModulationFactor
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.matplotlib_ import overlay_tag, save_current_figure


class TestModulationFactor(unittest.TestCase):

    """Unit test for xModulationFactor.
    """

    @classmethod
    def setUpClass(self):
        """Setup---here we essentially create the modulation factor.
        """
        self.measx, self.measy = numpy.loadtxt(GPD_MODF_FILE_PATH, unpack=True)
        file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.mrf')
        self.modf = xModulationFactor(file_path)

    def test_constant(self, num_events=1000000, interactive=False):
        """Test the modulation factor as a random number generator when
        both the polarization angle and degrees are energy- and
        time-independent.
        """
        polarization_angle = 20.
        polarization_degree = 1.
        self.modf.build_generator(polarization_angle, polarization_degree)
        self.modf.generator.plot(show=interactive)
        emin = self.modf.xmin()
        emax = self.modf.xmax()
        E = numpy.random.uniform(emin, emax, num_events)
        phi = self.modf.rvs(E)
        ebinning = numpy.linspace(emin, emax, 10)
        phi_binning = numpy.linspace(0, 360, 60)
        fit_results = []
        for _emin, _emax in zip(ebinning[:-1], ebinning[1:]):
            _emean = 0.5*(_emin + _emax)
            _mask = (E > _emin)*(E < _emax)
            _phi = phi[_mask]
            _hist = plt.hist(_phi, bins=phi_binning, histtype='step')
            _fr = xModulationFactor.fit_histogram(_hist)
            _fr.emean = _emean
            fit_results.append(_fr)
            _fr.plot()
            plt.axis([0., 360., 0., 1.2*_hist[0].max()])
            overlay_tag()
            plt.text(0.1, 0.1, 'Energy: %.2f--%.2f keV' % (_emin, _emax),
                     transform=plt.gca().transAxes)
            plt.text(0.1, 0.05, str(fit_results), transform=plt.gca().transAxes)
        _x = [_fr.emean for _fr in fit_results]
        _y = [_fr.phi0 for _fr in fit_results]
        _dy = [_fr.phi0_err for _fr in fit_results]
        plt.clf()
        plt.errorbar(_x, _y, yerr=_dy, fmt='o')
        plt.axis([emin, emax, polarization_angle-5, polarization_angle+5])
        plt.xlabel('Energy [keV]')
        plt.ylabel('Modulation angle [$^\circ$]')
        #plt.show()
        plt.clf()
        _y = [_fr.vis for _fr in fit_results]
        _dy = [_fr.vis_err for _fr in fit_results]
        plt.errorbar(_x, _y, yerr=_dy, fmt='o')
        plt.axis([emin, emax, 0, 1])
        self.modf.plot(show=interactive)
        plt.xlabel('Energy [keV]')
        plt.ylabel('Modulation visibility')
        #plt.show()

if __name__ == '__main__':
    unittest.main()
