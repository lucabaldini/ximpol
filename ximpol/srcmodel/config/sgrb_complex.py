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
from ximpol.srcmodel.roi import xUniformDisk, xROIModel
from ximpol.srcmodel.polarization import constant
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear


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
    fmt = dict(xname='Energy', xunits='keV', yname='Flux',
               yunits='cm$^{-2}$ s$^{-1}$ keV$^{-1}$')
    return xInterpolatedUnivariateSplineLinear(energy, flux, **fmt)


ROI_MODEL = xROIModel(266.8, -28.46)

spectral_model_b1 = parse_spectral_model('spec_model_SgrB1.txt')

def energy_spectrum_b1(E, t):
    """
    """
    return spectral_model_b1(E)

polarization_degree_b1 = constant(0.405)
polarization_angle_b1 = constant(numpy.radians(88.3))
SgrB1 = xUniformDisk('Sgr B1', 266.75833, -28.5325, angular_radius(6.),
                     energy_spectrum_b1, polarization_degree_b1,
                     polarization_angle_b1)

spectral_model_b2 = parse_spectral_model('spec_model_SgrB2.txt')

def energy_spectrum_b2(E, t):
    """
    """
    return spectral_model_b2(E)

polarization_degree_b2 = constant(0.455)
polarization_angle_b2 = constant(numpy.radians(84.4))
SgrB2 = xUniformDisk('Sgr B2', 266.835, -28.38528, angular_radius(5.),
                     energy_spectrum_b2, polarization_degree_b2,
                     polarization_angle_b2)

ROI_MODEL.add_sources(SgrB1, SgrB2)


if __name__ == '__main__':
    print(ROI_MODEL)
    from ximpol.utils.matplotlib_ import pyplot as plt
    plt.figure('Sgr complex')
    spectral_model_b1.plot(show=False, label='Sgr B1')
    spectral_model_b2.plot(show=False, label='Sgr B2')
    plt.yscale('log')
    plt.axis([1, 10, 1e-7, 1e-3])
    plt.legend(bbox_to_anchor=(0.35, 0.85))
    plt.show()
