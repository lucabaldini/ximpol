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

from ximpol.srcmodel.roi import xUniformDisk, xROIModel
from ximpol.srcmodel.spectrum import power_law

RA = 10.
DEC = 10.
RAD = 0.03
POL_MAX = 0.8

ROI_MODEL = xROIModel(RA, DEC)

def polarization_degree(E, t, ra, dec):
    distance = numpy.sqrt((ra-RA)**2 + (dec-DEC)**2)
    return POL_MAX*(distance/RAD)

def polarization_angle(E, t, ra, dec):
    return numpy.arctan2(ra-RA,dec-DEC)

energy_spectrum = power_law(10., 2.)
disk = xUniformDisk('disk', RA, DEC, RAD, energy_spectrum,
                    polarization_degree, polarization_angle)

ROI_MODEL.add_source(disk)


if __name__ == '__main__':
    print(ROI_MODEL)
