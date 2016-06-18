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
from ximpol.srcmodel.polarization import constant

POL_DEGREE_BKG = constant(0)
POL_ANGLE_BKG = constant(numpy.radians(0))

POL_DEGREE_JET = constant(0.2)
POL_ANGLE_JET = constant(numpy.radians(70))

POL_DEGREE_CORE = constant(0)
POL_ANGLE_CORE = constant(numpy.radians(0))

POLARIZATION_DICT = {0: [POL_DEGREE_BKG, POL_ANGLE_BKG],
                     1: [POL_DEGREE_JET, POL_ANGLE_JET],
                     2: [POL_DEGREE_CORE,POL_ANGLE_CORE]}

if __name__ == '__main__':
    print 'Configuration file M87 for Chandra-to-ximpol converter.'
