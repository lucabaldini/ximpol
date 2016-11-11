#!/usr/bin/env python
#
# Copyright (C) 2015--2016, the ximpol team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU GengReral Public Licensese as published by
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
import sys
import os

from ximpol.srcmodel.gabs import xpeInterstellarAbsorptionModel
from ximpol.srcmodel.roi import xPointSource, xROIModel
from ximpol.srcmodel.polarization import constant
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol import XIMPOL_CONFIG
from ximpol.utils.logging_ import logger
from ximpol.utils.units_ import keV2erg

GRS1915_RA = 288.7958
GRS1915_DEC = 10.9456

MIN_ENERGY = 0.5  #keV
MAX_ENERGY = 10.  #keV
FLUX = 1.0e-8     #erg/s/cm2
GAL_NH = 1.39e22  #cm-2
SPIN = 0.998      #0, 0.5, 0.9, 0.998

FILE_PATH = os.path.join(XIMPOL_CONFIG, 'ascii',
                         'bh_spin/grs1915_105_a%s_t1_i70.txt'%(SPIN))
print "Using file ",FILE_PATH

def scale_flux(energy, flux, emin=2., emax=8.):
    ism_model = xpeInterstellarAbsorptionModel()
    ism_trans = ism_model.transmission_factor(GAL_NH)
    flux *= ism_trans(energy)
    norm = xInterpolatedUnivariateSplineLinear(energy,energy*flux).integral(
                                                                    emin, emax)
    norm = keV2erg(norm)
    scale_factor =  FLUX/norm
    return flux*scale_factor

def polarization_degree(E, t, ra, dec):
    return pol_degree_spline(E)

def polarization_angle(E, t, ra, dec):
    return numpy.radians(pol_angle_spline(E))

def energy_spectrum(E, t):
    return _energy_spectrum(E)

#load all the relevant information from the txt file
_energy, _flux, _Q_E, _U_E, _V_E, _pol_degree, _pol_angle,_cpol_angle  = numpy.loadtxt(FILE_PATH, unpack=True)

_mask = (_energy >= MIN_ENERGY)*(_energy <= MAX_ENERGY)
_energy = _energy[_mask]
_pol_degree = _pol_degree[_mask]

fmt = dict(xname='Energy [keV]', yname='Polarization degree')
pol_degree_spline = xInterpolatedUnivariateSplineLinear(_energy, _pol_degree,
                                                        **fmt)

_pol_angle = 90. + _pol_angle[_mask]
fmt = dict(xname='Energy [keV]', yname='Polarization angle [$^\\circ$]')
pol_angle_spline = xInterpolatedUnivariateSplineLinear(_energy, _pol_angle,
                                                       **fmt)

flux = _flux[_mask]
_flux = scale_flux(_energy, flux)
fmt = dict(xname='Energy', xunits='keV', yname='Flux',
               yunits='cm$^{-2}$ s$^{-1}$ keV$^{-1}$')
_energy_spectrum = xInterpolatedUnivariateSplineLinear(_energy, _flux, **fmt)

ROI_MODEL = xROIModel(GRS1915_RA, GRS1915_DEC)

src = xPointSource('GRS1915_105', ROI_MODEL.ra, ROI_MODEL.dec, energy_spectrum,
                   polarization_degree, polarization_angle)

ROI_MODEL.add_source(src)

if __name__ == '__main__':
    print(ROI_MODEL)
    from ximpol.utils.matplotlib_ import pyplot as plt
    fig = plt.figure('Spectrum')
    _energy_spectrum.plot(show=False, logy=True)
    fig = plt.figure('Polarization angle')
    pol_angle_spline.plot(show=False)
    fig = plt.figure('Polarization degree')
    pol_degree_spline.plot(show=False)
    plt.show()
