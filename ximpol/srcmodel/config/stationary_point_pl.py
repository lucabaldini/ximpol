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

from ximpol.srcmodel.source import xPointSource


def dNdE(E, t):
    """Photon energy spectrum as a function of energy and time.
    """
    return 10.*numpy.power(E, -2.)

def polarization_degree(E, t):
    """Polarization degree as a function of energy and time.
    """
    return 1.

def polarization_angle(E, t):
    """Polarization angle as a function of energy and time.
    """
    return numpy.radians(44.)


source = xPointSource(name='test point source', ra=10., dec=10.)
source.spectrum = dNdE
source.polarization_degree = polarization_degree
source.polarization_angle = polarization_angle
source.identifier = 1
