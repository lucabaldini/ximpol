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


import os
import numpy

from ximpol import XIMPOL_DOC_FIGURES
from ximpol.irf import load_irfs
from ximpol.utils.matplotlib_ import save_current_figure
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.core.spline import xInterpolatedBivariateSplineLinear
from ximpol.srcmodel.spectrum import power_law, xCountSpectrum


IRF_NAME = 'xipe_baseline'
T_MIN = 0.
T_MAX = 1000.
OUTPUT_FOLDER = XIMPOL_DOC_FIGURES

aeff, psf, modf, edisp = load_irfs(IRF_NAME)

def C(t):
    return 1. - 0.9*t/(T_MAX - T_MIN)

def Gamma(t):
    return 1.5 + t/(T_MAX - T_MIN)

spectrum = power_law(C, Gamma)

plt.figure('Source spectrum')
_t = numpy.linspace(T_MIN, T_MAX, 100)
_e = aeff.x
fmt = dict(xname='Time', xunits='s', yname='Energy', yunits='keV',
           zname='dN/dE', zunits='cm$^{-2}$ s$^{-1}$ keV$^{-1}$')
spectrum_spline = xInterpolatedBivariateSplineLinear(_t, _e, spectrum, **fmt)
spectrum_spline.plot(show=False, logz=True)
save_current_figure('source_spectrum.png', OUTPUT_FOLDER, clear=False)

plt.figure('Source spectrum slices')
for t in [T_MIN, 0.5*(T_MIN + T_MAX), T_MAX]:
    spectrum_spline.vslice(t).plot(show=False, logx=True, logy=True,
                                   label='t = %d s' % t)
plt.legend(bbox_to_anchor=(0.75, 0.95))
plt.axis([1, 10, None, None])
save_current_figure('source_spectrum_slices.png', OUTPUT_FOLDER, clear=False)

plt.figure('Count spectrum')
count_spectrum = xCountSpectrum(spectrum, aeff, _t)
count_spectrum.plot(show=False, logz=True)
save_current_figure('count_spectrum.png', OUTPUT_FOLDER, clear=False)

plt.figure('Count spectrum slices')
for t in [T_MIN, 0.5*(T_MIN + T_MAX), T_MAX]:
    count_spectrum.vslice(t).plot(show=False, logx=True, logy=True,
                                  label='t = %d s' % t)
plt.legend(bbox_to_anchor=(0.9, 0.95))
plt.axis([1, 10, None, None])
save_current_figure('count_spectrum_slices.png', OUTPUT_FOLDER, clear=False)

plt.figure('Light curve')
count_spectrum.light_curve.plot(show=False)
save_current_figure('light_curve.png', OUTPUT_FOLDER, clear=False)

plt.figure('Light curve times')
count_spectrum.light_curve.plot(show=False)
for t in count_spectrum.light_curve.rvs(100):
    plt.plot([t, t], [0, 20], color='gray', linestyle='-', linewidth=2)
save_current_figure('light_curve_random.png', OUTPUT_FOLDER, clear=False)

plt.show()
