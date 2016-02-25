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
from ximpol.srcmodel.spectrum import power_law, constant


ROI_MODEL = xROIModel(10., 10.)

energy_spectrum1 = power_law(10., 2.)
polarization_degree1 = constant(1.0)
polarization_angle1 = constant(numpy.radians(65.))
src1 = xPointSource('Point source 1', 9.9750, 9.9833, energy_spectrum1,
                    polarization_degree1, polarization_angle1)

energy_spectrum2 = power_law(15., 3.)
polarization_degree2 = constant(0.0)
polarization_angle2 = constant(numpy.radians(0))
src2 = xPointSource('Point source 2', 10., 10., energy_spectrum2,
                    polarization_degree2, polarization_angle2)

ROI_MODEL.add_sources(src1, src2)


if __name__ == '__main__':
    print(ROI_MODEL)
