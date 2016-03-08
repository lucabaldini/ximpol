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
from ximpol.srcmodel.polarization import xPolMap, constant
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.utils.logging_ import logger
from ximpol import XIMPOL_CONFIG


ROI_MODEL = xROIModel(350.8664167, 58.8117778)

total_spec_file_path = os.path.join(XIMPOL_CONFIG, 'ascii',
                                    'CasA_total_spectrum.txt')
nonthermal_spec_file_path = os.path.join(XIMPOL_CONFIG, 'ascii',
                                         'CasA_nonthermal_spectrum.txt')

le_img_file_path = os.path.join(XIMPOL_CONFIG, 'fits', 'casa_1p5_3p0_keV.fits')
he_img_file_path = os.path.join(XIMPOL_CONFIG, 'fits', 'casa_4p0_6p0_keV.fits')

def parse_spectrum(file_path, emin=1., emax=12.):
    """Parse the input file with the spectral point.
    """
    logger.info('Parsing input file %s...' % file_path)
    energy, flux = numpy.loadtxt(file_path, unpack=True)
    _mask = (energy >= emin)*(energy <= emax)
    return energy[_mask], flux[_mask]

# Parse the total spectrum
_energy, _flux = parse_spectrum(total_spec_file_path)
fmt = dict(xname='Energy', xunits='keV', yname='Flux',
           yunits='cm$^{-2}$ s$^{-1}$ keV$^{-1}$')
total_spectral_model = xInterpolatedUnivariateSplineLinear(_energy, _flux,
                                                           **fmt)
# Parse the non-thermal spectrum
_energy, _flux = parse_spectrum(nonthermal_spec_file_path)
fmt = dict(xname='Energy', xunits='keV', yname='Flux',
           yunits='cm$^{-2}$ s$^{-1}$ keV$^{-1}$')
nonthermal_spectral_model = xInterpolatedUnivariateSplineLinear(_energy, _flux,
                                                                **fmt)
# Subtract the two to get the thermal component.
thermal_spectral_model = total_spectral_model - nonthermal_spectral_model

def thermal_energy_spectrum(E, t):
    return thermal_spectral_model(E)

def nonthermal_energy_spectrum(E, t):
    return nonthermal_spectral_model(E)

thermal_polarization_angle = constant(0.)
thermal_polarization_degree = constant(0.)

# Read the polarization maps for the non-thermal component.
pol_mapx_path = os.path.join(XIMPOL_CONFIG, 'fits', 'casa_pol_x.fits')
pol_mapy_path = os.path.join(XIMPOL_CONFIG, 'fits', 'casa_pol_y.fits')
polarization_map = xPolMap(pol_mapx_path, pol_mapy_path)

def nonthermal_polarization_angle(E, t, ra, dec):
    return polarization_map.polarization_angle(ra, dec)

def nonthermal_polarization_degree(E, t, ra, dec):
    return polarization_map.polarization_degree(ra, dec)

thermal_component = xExtendedSource('Cas A thermal', le_img_file_path,
                                    thermal_energy_spectrum,
                                    thermal_polarization_degree,
                                    thermal_polarization_angle)
nonthermal_component = xExtendedSource('Cas A non-thermal', he_img_file_path,
                                       nonthermal_energy_spectrum,
                                       nonthermal_polarization_degree,
                                       nonthermal_polarization_angle)

ROI_MODEL.add_sources(thermal_component, nonthermal_component)


if __name__ == '__main__':
    from ximpol.utils.matplotlib_ import pyplot as plt
    print(ROI_MODEL)
    fig = plt.figure('Energy spectrum')
    total_spectral_model.plot(logy=True, show=False)
    nonthermal_spectral_model.plot(logy=True, show=False)
    thermal_spectral_model.plot(logy=True, show=False)
    plt.show()
