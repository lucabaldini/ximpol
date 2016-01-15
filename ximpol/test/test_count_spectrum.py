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
from ximpol.irf.arf import xEffectiveArea
from ximpol.srcmodel.spectrum import xCountSpectrum
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.matplotlib_ import overlay_tag, save_current_figure


class TestCountSpectrum(unittest.TestCase):

    """Unit test for xCountSpectrum.
    """

    @classmethod
    def setUpClass(self):
        """Setup---here we essentially create the effective area.

        >>> file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.arf')
        >>> self.aeff = xEffectiveArea(file_path)
        """
        file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.arf')
        self.aeff = xEffectiveArea(file_path)

    def test_power_law_stationary(self):
        """Test a time-independent power law.

        This creates a count spectrum with no time dependence, i.e., with
        only two (identical) interpolating time points on the auxiliary
        (time) axis. The power-law parameters (C and gamma) are constant

        >>> tmin = 0.
        >>> tmax = 100.
        >>> C = 1.
        >>> Gamma = 2.
        >>>
        >>> def powerlaw(E, t):
        >>>     return C*numpy.power(E, -Gamma)
        >>>
        >>> _t = numpy.linspace(tmin, tmax, 2)
        >>> count_spectrum = xCountSpectrum(powerlaw, self.aeff, _t)

        and the underlying xUnivariateAuxGenerator looks like.

        .. image:: ../figures/test_power_law_stationary_2d.png

        Then a vertical slice (i.e., an interpolated linear spline) is taken
        in the middle of the auxiliary axis

        >>> tref = 0.5*(tmin + tmax)
        >>> ref_slice = count_spectrum.slice(tref)

        and the y-values of the spline are compared with the direct product of
        the effective area and the input spectrum:

        >>> _x = self.aeff.x
        >>> _y = C*numpy.power(_x, -Gamma)*self.aeff.y

        (Note that in general the two are technically not the same thing, as
        going from the count spectrum to the slice we do interpolate in time,
        although in this particular case the interpolation is trivial).
        If everything goes well, they should be on top of each other.
        The figure below is also showing the original power-law spectrum
        multiplied by the peak effective area.

        .. image:: ../figures/test_power_law_stationary_slice.png

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
        count_spectrum.plot(show=False)
        overlay_tag(color='white')
        save_current_figure('test_power_law_stationary_2d.png')

        ref_slice = count_spectrum.slice(tref)
        _x = self.aeff.x
        _y = C*numpy.power(_x, -Gamma)*self.aeff.y.max()
        plt.plot(_x, _y, '-', label='Original power-law spectrum')
        _y = C*numpy.power(_x, -Gamma)*self.aeff.y
        delta = abs((_y - ref_slice(_x))/_y).max()
        #self.assertTrue(delta < 1e-3, 'max deviation %.9f' % delta)
        plt.plot(_x, _y, 'o', label='Direct convolution with aeff')
        ref_slice.plot(logx=True, logy=True, show=False,
                       label='xCountSpectrum output')
        overlay_tag()
        plt.text(0.1, 0.1, 'Max. difference = %.3e' % delta,
                 transform=plt.gca().transAxes)
        plt.legend(bbox_to_anchor=(0.65, 0.5))
        save_current_figure('test_power_law_stationary_slice.png')

    def test_power_law_variable(self):
        """Test a time-dependent power law.

        This creates a time-dependent count spectrum, where the two parameters
        of the underlying power law (C and Gamma) vary linearly with time, in
        opposite direction, between 1 and 2.

        >>> def C(t):
        >>>     return 1. + (t - tmin)/(tmax - tmin)
        >>>
        >>> def Gamma(t):
        >>>    return 2. - (t - tmin)/(tmax - tmin)
        >>>
        >>> def powerlaw(E, t):
        >>>    return C(t)*numpy.power(E, -Gamma(t))

        (Beware: this does not mean that you can interpolate linearly between
        the two time extremes, as both parameters vary at the same time and
        the spectral shape does not evolve linearly with time---we're sampling
        the time axis with 100 points).
        The underlying xUnivariateAuxGenerator looks like.

        .. image:: ../figures/test_power_law_variable_2d.png

        Then a vertical slice (i.e., an interpolated linear spline) is taken
        in the middle of the auxiliary axis and the y-values of the spline
        are compared with the direct product of the effective area and the
        count spectrum (evaluated at the same time). If everything goes well,
        they should be on top of each other. The figure below is also showing
        the orignal power-law spectrum multiplied by the peak effective area.

        .. image:: ../figures/test_power_law_variable_slice.png

        Finally, we do test the light-curve building by comparing it with the
        values from a direct intergration of the vertical slices on a
        fixed-spacing grid. Note that, since the normalization increases with
        time and the spectral index becomes harder, the light-curve increases
        more than linearly.

        .. image:: ../figures/test_power_law_variable_lc.png

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
        count_spectrum.plot(show=False)
        overlay_tag(color='white')
        save_current_figure('test_power_law_variable_2d.png')

        ref_slice = count_spectrum.slice(tref)
        _x = self.aeff.x
        _y = self.aeff.y.max()*C(tref)*numpy.power(_x, -Gamma(tref))
        plt.plot(_x, _y, '-', label='Original power-law spectrum')
        _y = C(tref)*numpy.power(_x, -Gamma(tref))*self.aeff.y
        delta = abs((_y - ref_slice(_x))/_y).max()
        #self.assertTrue(delta < 1e-3, 'max deviation %.9f' % delta)
        plt.plot(_x, _y, 'o', label='Direct convolution with aeff')
        ref_slice.plot(logx=True, logy=True, show=False,
                       label='xCountSpectrum output')
        overlay_tag()
        plt.text(0.1, 0.1, 'Max. difference = %.3e' % delta,
                 transform=plt.gca().transAxes)
        plt.legend(bbox_to_anchor=(0.65, 0.5))
        save_current_figure('test_power_law_variable_slice.png')

        _x = numpy.linspace(tmin, tmax, 33)
        _y = []
        for _xp in _x:
            _y.append(count_spectrum.slice(_xp).norm())
        _y = numpy.array(_y)
        plt.plot(_x, _y, 'o', label='Direct integral flux values')
        delta = abs((_y - count_spectrum.light_curve(_x))/_y).max()
        self.assertTrue(delta < 1e-3, 'max deviation %.9f' % delta)
        count_spectrum.light_curve.plot(show=False,
                                        label='xCountSpectrum light-curve')
        overlay_tag()
        plt.legend(bbox_to_anchor=(0.65, 0.75))
        plt.text(0.5, 0.1, 'Max. difference = %.3e' % delta,
                 transform=plt.gca().transAxes)
        save_current_figure('test_power_law_variable_lc.png')

    def test_power_law_rvs(self, num_events=1000000):
        """Test the generation of event energies from a count power-law
        spectrum convoluted with the effective area.

        This turned out to be more tricky than we anticipated. Since the
        convolution of the source spectrum with the effective area falls
        pretty quickly at high energy, there's typically very few events above
        a few keV and, as a consequence, the slope of the corresponding ppf is
        fairly steep close to one.

        .. image:: ../figures/test_power_law_rvs_vppf.png

        This implies that the ppf in each slice must be properly sampled
        close to 1 (initial tests showed, e.g., that a uniform grid with
        100 points between 0 and 1 was not enough to throw meaningful random
        numbers for a typical power-law source spectrum). This is particularly
        true for soft spectral indices---which is why we picked `Gamma = 3.`
        for this test.

        .. image:: ../figures/test_power_law_rvs_counts.png

        """
        tmin = 0.
        tmax = 100.
        tref = 0.5*(tmin + tmax)
        C = 1.
        Gamma = 3.

        def powerlaw(E, t):
            """Function defining a time-dependent energy spectrum.
            """
            return C*numpy.power(E, -Gamma)

        _t = numpy.linspace(tmin, tmax, 2)
        count_spectrum = xCountSpectrum(powerlaw, self.aeff, _t)

        count_spectrum.vppf.vslice(tref).plot(show=False, overlay=True)
        overlay_tag()
        save_current_figure('test_power_law_rvs_vppf.png')

        ref_slice = count_spectrum.slice(tref)
        _time = numpy.zeros(num_events)
        _time.fill(tref)
        _energy = count_spectrum.rvs(_time)
        _binning = numpy.linspace(self.aeff.xmin(), self.aeff.xmax(), 100)
        obs, bins, patches = plt.hist(_energy, bins=_binning,
                                      histtype='step', label='Random energies')
        plt.yscale('log')
        # We want to overlay the reference count-spectrum slice, normalized
        # to the total number of events simulated.
        bin_width = (bins[1] - bins[0])
        scale = num_events*bin_width/ref_slice.norm()
        _x = 0.5*(bins[:-1] + bins[1:])
        _y = scale*ref_slice(_x)
        plt.plot(_x, _y, label='Underlying pdf')
        # And, for the chisquare, we do correctly integrate the slice in each
        # energy bin, rather than evaluating it at the bin center.
        exp = []
        scale = num_events/ref_slice.norm()
        for _emin, _emax in zip(bins[:-1], bins[1:]):
            exp.append(scale*ref_slice.integral(_emin, _emax))
        exp = numpy.array(exp)
        chisquare = ((exp - obs)**2/exp).sum()
        plt.text(0.5, 0.1, '$\chi^2$/ndof = %.2f/%d' % (chisquare, len(obs)),
                 transform=plt.gca().transAxes)
        plt.legend(bbox_to_anchor=(0.5, 0.5))
        overlay_tag()
        save_current_figure('test_power_law_rvs_counts.png')


if __name__ == '__main__':
    unittest.main()
