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

from ximpol.srcmodel.roi import xPointSource, xROIModel
from ximpol.srcmodel.spectrum import power_law, constant


ROI_MODEL = xROIModel(10., 10.)

source1 = xPointSource(name='source1', ra=9.9750, dec=9.9833)
source1.spectrum = power_law(10., 2.)
source1.polarization_degree = constant(1.0)
source1.polarization_angle = constant(numpy.radians(65.))

source2 = xPointSource(name='source2', ra=10., dec=10.)
source2.spectrum = power_law(15., 3.)
source2.polarization_degree = constant(0.0)
source2.polarization_angle = constant(numpy.radians(0))

ROI_MODEL.add_sources(source1, source2)


if __name__ == '__main__':
    print(ROI_MODEL)
