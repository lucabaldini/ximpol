#!/usr/bin/env python
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


import numpy
import os

from ximpol.srcmodel.roi import xExtendedSource, xROIModel
from ximpol.srcmodel.spectrum import power_law
from ximpol.srcmodel.polarization import constant
from ximpol import XIMPOL_CONFIG


ROI_MODEL = xROIModel(83.633083, 22.014500)

img_file_path = os.path.join(XIMPOL_CONFIG, 'fits', 'crab_0p3_10p0_keV.fits')
energy_spectrum = power_law(9.59, 2.108)
polarization_degree = constant(0.157)
polarization_angle = constant(numpy.radians(161.1))
crab_nebula = xExtendedSource('Crab nebula', img_file_path, energy_spectrum,
                              polarization_degree, polarization_angle)

ROI_MODEL.add_source(crab_nebula)


if __name__ == '__main__':
    print(ROI_MODEL)
