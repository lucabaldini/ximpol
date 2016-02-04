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

from ximpol import XIMPOL_SRCMODEL
from ximpol.srcmodel.roi import xUniformDisk, xPointSource, xROIModel
from ximpol.srcmodel.spectrum import constant, xTabulatedStationarySpectrum


def angular_radius(physical_radius):
    """Return the angular radius (in degrees) subtended by a cloud located at
    the Galactic Center (assumed to be 8 kpc away from the Earth), given the
    physical extension radius in pc.
    """
    return numpy.degrees(numpy.arctan(physical_radius/8000.))

def parse_spectral_model(file_name):
    """Parse a spectral model written by XSPEC and return the corresponding
    interpolated univariate spline.
    """
    file_path = os.path.join(XIMPOL_SRCMODEL, 'ascii', file_name)
    data = numpy.loadtxt(file_path, unpack=True)
    energy, flux = data[0], data[2]
    return xTabulatedStationarySpectrum(energy, flux)


ROI_MODEL = xROIModel(266.8, -28.44)

SgrB1 = xUniformDisk('Sgr B1', 266.75833, -28.5325, angular_radius(6.))
SgrB1.spectrum = parse_spectral_model('spec_model_SgrB1.txt')
SgrB1.polarization_degree = constant(0.405)
SgrB1.polarization_angle = constant(numpy.radians(88.3))

SgrB2 = xUniformDisk('Sgr B2', 266.835, -28.38528, angular_radius(5.))
SgrB2.spectrum = parse_spectral_model('spec_model_SgrB2.txt')
SgrB2.polarization_degree = constant(0.455)
SgrB2.polarization_angle = constant(numpy.radians(84.4))

ROI_MODEL.add_sources(SgrB1, SgrB2)
