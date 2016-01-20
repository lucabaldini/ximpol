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
import sys

from ximpol.irf.mrf import xAzimuthalResponseGenerator
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.matplotlib_ import overlay_tag, save_current_figure


"""We explictely set the random seed to have reproducible results.
"""
numpy.random.seed(0)


class TestAzimuthalResponse(unittest.TestCase):

    """Unit test for xAzimuthalResponseGenerator
    """

    @classmethod
    def setUpClass(self):
        """Setup---here we essentially create the effective area.
        """
        self.generator = xAzimuthalResponseGenerator()
        self.interactive = sys.flags.interactive

    def test_pdf(self):
        """Test the one-dimensional azimuthal response underlying pdf.
        """
        self.generator.plot(show=self.interactive)
        overlay_tag(color='white')
        save_current_figure('test_azimuthal_resp_generator.png',
                            show=self.interactive)
        phi = numpy.linspace(0., 2*numpy.pi, 100)
        for visibility in numpy.linspace(1, 0, 5):
            pdf = self.generator.pdf(phi, visibility)
            plt.plot(phi, pdf, label='$\\xi = %.2f$' % visibility)
            spline = xInterpolatedUnivariateSplineLinear(phi, pdf)
            norm = spline.norm()
            self.assertTrue(abs(norm - 1.) < 1e-5,
                            'Normalization is %.3e' % norm)
        plt.axis([0., 2*numpy.pi, None, None])
        plt.xlabel('$\\phi$ [rad]')
        plt.ylabel('pdf($\\phi$) [1/rad]')
        plt.legend(bbox_to_anchor=(0.88, 0.92))
        overlay_tag()
        save_current_figure('test_azimuthal_resp_pdf.png',
                            show=self.interactive)

    def test_cdf(self):
        """Test the one-dimensional azimuthal response underlying pdf.
        """
        phi = numpy.linspace(0., 2*numpy.pi, 100)
        for visibility in numpy.linspace(1, 0, 5):
            cdf = self.generator.cdf(phi, visibility)
            plt.plot(phi, cdf, label='$\\xi = %.2f$' % visibility)
            spline = xInterpolatedUnivariateSplineLinear(phi, cdf)
            self.assertTrue(abs(spline(0.)) < 1e-5, 'cdf(0) = %.3e' % spline(0))
            self.assertTrue(abs(spline(2*numpy.pi) - 1) < 1e-5,
                            'cdf(2pi) = %.3e' % spline(2*numpy.pi))
        plt.axis([0., 2*numpy.pi, None, None])
        plt.xlabel('$\\phi$ [rad]')
        plt.ylabel('cdf($\\phi$)')
        plt.legend(bbox_to_anchor=(0.4, 0.92))
        overlay_tag()
        save_current_figure('test_azimuthal_resp_cdf.png',
                            show=self.interactive)

    def test_rvs(self):
        """Test the random number generation.
        """
        visibility = numpy.zeros(1000000)
        visibility.fill(0.5)
        phi = self.generator.rvs_phi(visibility, 0.25*numpy.pi)
        hist = plt.hist(phi, bins=numpy.linspace(0, 2*numpy.pi, 100),
                        histtype='step', label='Random angles')
        fit_results = self.generator.fit_histogram(hist)
        fit_results.plot()
        plt.xlabel('$\\phi$ [rad]')
        plt.axis([0, 2*numpy.pi, 0, None])
        overlay_tag()
        save_current_figure('test_azimuthal_resp_rvs.png',
                            show=self.interactive)


if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)
