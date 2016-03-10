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


"""Configuration file for a semi-realistic model of Cas A.

The spectral model is taken from
E.A. Helder and J. Vink, "Characterizing the non-thermal emission of Cas A",
Astrophys.J. 686 (2008) 1094--1102, http://arxiv.org/abs/0806.3748,
which seems to be one of the few instances where an actual spectrum in physical
units (i.e., not a count spectrum) is presented.
We grabbed by hand the data points and we're calling "thermal" whatever is in
the lines and "non-thermal" the rest.

We have two images of Cas A, at low (1.5--3.0 keV) and high (4.0--6.0 keV) and,
due to the absence of lines between 4 and 6 keV we're attaching the latter to
the non-thermal spectrum and the former to the thermal component.

The polarization map is a simple geometrical, radially-symmetric, model.
"""


def parse_spectral_model(file_name, emin=1., emax=15.):
    """Parse the input file with the spectral point.
    """
    file_path = os.path.join(XIMPOL_CONFIG, 'ascii', file_name)
    logger.info('Parsing input file %s...' % file_path)
    _energy, _flux = numpy.loadtxt(file_path, delimiter=',', unpack=True)
    _mask = (_energy >= emin)*(_energy <= emax)
    _energy = _energy[_mask]
    _flux = _flux[_mask]
    fmt = dict(xname='Energy', xunits='keV', yname='Flux',
               yunits='cm$^{-2}$ s$^{-1}$ keV$^{-1}$')
    return xInterpolatedUnivariateSplineLinear(_energy, _flux, **fmt)


ROI_MODEL = xROIModel(350.8664167, 58.8117778)

le_img_file_path = os.path.join(XIMPOL_CONFIG, 'fits', 'casa_1p5_3p0_keV.fits')
he_img_file_path = os.path.join(XIMPOL_CONFIG, 'fits', 'casa_4p0_6p0_keV.fits')

# Read in the spectral models.
total_spectral_model = parse_spectral_model('casa_total_spectrum.csv')
nonthermal_spectral_model = parse_spectral_model('casa_nonthermal_spectrum.csv')
thermal_spectral_model = total_spectral_model - nonthermal_spectral_model

def thermal_energy_spectrum(E, t):
    return thermal_spectral_model(E)

def nonthermal_energy_spectrum(E, t):
    return nonthermal_spectral_model(E)

# The thermal component is totally unpolarized.
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


def display():
    """
    """
    from ximpol.utils.matplotlib_ import pyplot as plt
    print(ROI_MODEL)
    fig = plt.figure('Energy spectrum')
    total_spectral_model.plot(logy=True, show=False, label='Total')
    nonthermal_spectral_model.plot(logy=True, show=False, label='Non-thermal')
    thermal_spectral_model.plot(logy=True, show=False, label='Thermal')
    plt.legend(bbox_to_anchor=(0.95, 0.95))
    fig = thermal_component.image.plot(show=False)
    fig.add_label(0.1, 0.92, '1.5-3 keV', relative=True, size='xx-large',
                  color='white', horizontalalignment='left')
    fig = nonthermal_component.image.plot(show=False)
    fig.add_label(0.1, 0.92, '4-6 keV', relative=True, size='xx-large',
                  color='white', horizontalalignment='left')
    plt.show()


if __name__ == '__main__':
    display()
