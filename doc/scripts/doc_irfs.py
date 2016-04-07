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

from ximpol import XIMPOL_DOC_FIGURES
from ximpol.irf import load_irfs
from ximpol.utils.matplotlib_ import save_current_figure
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear


IRF_NAME = 'xipe_baseline'
OUTPUT_FOLDER = XIMPOL_DOC_FIGURES

aeff, psf, modf, edisp = load_irfs(IRF_NAME)

plt.figure('On axis effective area')
xInterpolatedUnivariateSplineLinear.plot(aeff, show=False)
save_current_figure('aeff_on_axis.png', OUTPUT_FOLDER, clear=False)

plt.figure('Vignetting')
aeff.vignetting.plot(show=False)
save_current_figure('aeff_vignetting.png', OUTPUT_FOLDER, clear=False)

plt.figure('Edisp matrix')
edisp.matrix.plot(show=False)
save_current_figure('edisp_matrix.png', OUTPUT_FOLDER, clear=False)

plt.figure('Edisp slice')
_e = 6.
edisp.matrix.vslice(_e).plot(show=False, label='E = %.2f keV' % _e)
plt.axis([0, 256, None, None])
plt.legend(bbox_to_anchor=(0.45, 0.75))
save_current_figure('edisp_slice.png', OUTPUT_FOLDER, clear=False)

plt.figure('PSF')
psf.plot(show=False)
save_current_figure('psf.png', OUTPUT_FOLDER, clear=False)

plt.figure('Modulation factor')
modf.plot(show=False)
save_current_figure('modf.png', OUTPUT_FOLDER, clear=False)

plt.show()
