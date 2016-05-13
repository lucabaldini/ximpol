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

from ximpol.srcmodel.roi import xPointSource, xROIModel
from ximpol.srcmodel.spectrum import power_law
from ximpol.srcmodel.polarization import xPolarizationMap, constant
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.utils.logging_ import logger
from ximpol import XIMPOL_CONFIG


"""Configuration file for a model of Abell 478.
"""


def parse_spectral_model(file_name, emin=0.9, emax=11.):
    """Parse the input file with the spectral point.
    """
    file_path = os.path.join(XIMPOL_CONFIG, 'ascii', file_name)
    logger.info('Parsing input file %s...' % file_path)
    _energy, _flux, _fluxerr = numpy.loadtxt(file_path, unpack=True)
    _mask = (_energy >= emin)*(_energy <= emax)
    _energy = _energy[_mask]
    _flux = _flux[_mask]
    fmt = dict(xname='Energy', xunits='keV', yname='Flux',
               yunits='cm$^{-2}$ s$^{-1}$ keV$^{-1}$')
    return xInterpolatedUnivariateSplineLinear(_energy, _flux, **fmt)


ROI_MODEL = xROIModel(63.3363, 10.4764)

# Read in the spectral models.
spectral_model_spline = parse_spectral_model('Abell_478.txt')

def energy_spectrum(E, t):
    return spectral_model_spline(E)

polarization_degree = constant(0.)
polarization_angle = constant(0.)

abell478 = xPointSource('Abell 478', ROI_MODEL.ra, ROI_MODEL.dec,
                       energy_spectrum, polarization_degree, polarization_angle)

ROI_MODEL.add_source(abell478)


def display():
    """Display the source model.
    """
    from ximpol.utils.matplotlib_ import pyplot as plt
    from ximpol.srcmodel.img import xFITSImage

    print(ROI_MODEL)
    fig = plt.figure('Energy spectrum')
    spectral_model_spline.plot(logy=True, show=False, label='Total')
    plt.show()


if __name__ == '__main__':
    display()
