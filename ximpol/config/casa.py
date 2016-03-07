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
from ximpol.srcmodel.polarization import xPolMap
from ximpol import XIMPOL_CONFIG


ROI_MODEL = xROIModel(350.8664167, 58.8117778)

img_file_path = os.path.join(XIMPOL_CONFIG, 'fits', 'casa_1p5_3p0_keV.fits')

# This is actually for the Crab and needs to be fixed.
energy_spectrum = power_law(10., 2.)

# Read the polarization maps/
pol_mapx_path = os.path.join(XIMPOL_CONFIG, 'fits', 'casa_pol_x.fits')
pol_mapy_path = os.path.join(XIMPOL_CONFIG, 'fits', 'casa_pol_y.fits')
polarization_map = xPolMap(pol_mapx_path, pol_mapy_path)

def polarization_angle(E, t, ra, dec):
    return polarization_map.polarization_angle(ra, dec)

def polarization_degree(E, t, ra, dec):
    return polarization_map.polarization_degree(ra, dec)

casa = xExtendedSource('Cas A', img_file_path, energy_spectrum,
                       polarization_degree, polarization_angle)

ROI_MODEL.add_source(casa)


if __name__ == '__main__':
    print(ROI_MODEL)
    north = 350.84852, 58.839822
    west = 350.91201, 58.804452
    center = 350.85707, 58.81634
    for ra, dec in [north, west, center]:
        print ra, dec,\
            polarization_map.polarization_degree(ra, dec),\
            numpy.degrees(polarization_map.polarization_angle(ra, dec))
