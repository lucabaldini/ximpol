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


import numpy
import os

from ximpol.srcmodel.source import xExtendedSource, xROIModel
from ximpol.srcmodel.spectrum import power_law, constant

from ximpol import XIMPOL_SRCMODEL

ROI_MODEL = xROIModel(10., 10., 1.0)

img_file_path = os.path.join(XIMPOL_SRCMODEL, 'fits', 'crab_0p3_10p0_keV.fits')
source = xExtendedSource('Crab Nebula', img_file_path)
source.spectrum = power_law(10., 2.)
source.polarization_degree = constant(0.157)
source.polarization_angle = constant(numpy.radians(161.1))

ROI_MODEL.add_source(source)