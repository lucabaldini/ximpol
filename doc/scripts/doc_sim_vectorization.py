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
from ximpol.core.rand import xUnivariateGenerator
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.core.spline import xInterpolatedBivariateSplineLinear
from ximpol.srcmodel.spectrum import power_law, xCountSpectrum


OUTPUT_FOLDER = XIMPOL_DOC_FIGURES

_x = numpy.linspace(-5., 5., 100)
_y = (1./numpy.sqrt(2*numpy.pi))*numpy.exp(-_x**2/2.)
pdf = xUnivariateGenerator(_x, _y)
print pdf.norm()

plt.figure('pdf')
pdf.plot(show=False)

plt.figure('ppf')
pdf.ppf.plot(show=False)

print numpy.random.sample(10)

plt.show()
