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


"""Configuration file for a model of Abell 2142.
"""


def parse_spectral_model(file_name, emin=0.5, emax=15.):
    """Parse the input file with the spectral point.
    """
    file_path = os.path.join(XIMPOL_CONFIG, 'ascii', file_name)
    logger.info('Parsing input file %s...' % file_path)
    _energy, _binw, _flux, _fluxerr, _mod = numpy.loadtxt(file_path, unpack=True)
    _mask = (_energy >= emin)*(_energy <= emax)
    _energy = _energy[_mask]
    _mod = _mod[_mask]
    _mod /= _energy**2.
    fmt = dict(xname='Energy', xunits='keV', yname='Flux',
               yunits='cm$^{-2}$ s$^{-1}$ keV$^{-1}$')
    return xInterpolatedUnivariateSplineLinear(_energy, _mod, **fmt)


ROI_MODEL = xROIModel(239.5833, 27.2334)

# Read in the spectral models.
spectral_model_spline = parse_spectral_model('Abell_2142.txt')

def energy_spectrum(E, t):
    return spectral_model_spline(E)

polarization_degree = constant(0.)
polarization_angle = constant(0.)

abell2142 = xPointSource('Abell 2142', ROI_MODEL.ra, ROI_MODEL.dec,
                       energy_spectrum, polarization_degree, polarization_angle)

ROI_MODEL.add_source(abell2142)


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
