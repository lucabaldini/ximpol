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

from ximpol.srcmodel.roi import xPointSource, xROIModel
from ximpol.srcmodel.spectrum import power_law
from ximpol.srcmodel.polarization import constant


ROI_MODEL = xROIModel(90., 0.)
PL_NORM = 10.
PL_INDEX = 2.0
POL_DEGREE = 0.5
POL_ANGLE = numpy.radians(65.)

energy_spectrum = power_law(PL_NORM, PL_INDEX)
polarization_degree = constant(POL_DEGREE)
polarization_angle = constant(POL_ANGLE)

src = xPointSource('Point source', ROI_MODEL.ra, ROI_MODEL.dec, energy_spectrum,
                   polarization_degree, polarization_angle)

ROI_MODEL.add_source(src)


if __name__ == '__main__':
    print(ROI_MODEL)
