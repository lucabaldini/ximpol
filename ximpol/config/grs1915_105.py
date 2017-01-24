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
import scipy
import scipy.signal

from ximpol.srcmodel.roi import xPointSource, xROIModel
from ximpol.srcmodel.spectrum import power_law
from ximpol.srcmodel.polarization import constant
from ximpol.core.spline import xInterpolatedUnivariateSpline
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol import XIMPOL_CONFIG
from ximpol.utils.logging_ import logger


GRS1915_RA = 288.7958
GRS1915_DEC = 10.9456 

MIN_ENERGY = 0.1
MAX_ENERGY = 15.

spindegree = 0.998

#66 degree inclination
POL_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'ascii', 'bh_spin/Polarization_spin%s.txt'%(spindegree))

print "Using path",POL_FILE_PATH
FLUX_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'ascii', 'bh_spin/Flux_spin%s.txt'%(spindegree))


def polarization_degree(E, t, ra, dec):
    return pol_degree_spline(E)

def polarization_angle(E, t, ra, dec):
    return pol_angle_spline(E)

def energy_spectrum(E, t):
    return _energy_spectrum(E)


# Build the polarization degree as a function of the energy.
_energy, _pol_degree, _pol_angle = numpy.loadtxt(POL_FILE_PATH, unpack=True)

#Switch to have degrees
_pol_degree /= 100.
#_pol_degree = _pol_degree

_mask = (_energy >= MIN_ENERGY)*(_energy <= MAX_ENERGY)
_energy = _energy[_mask]
_pol_degree = _pol_degree[_mask]

fmt = dict(xname='Energy', yname='Polarization degree')
pol_degree_spline = xInterpolatedUnivariateSpline(_energy, _pol_degree, k=1, **fmt)
print "Pol angle in degrees",_pol_angle
#Pol angle in radians
_pol_angle = numpy.deg2rad(_pol_angle)

#Switched to have degrees and not radians
#_pol_angle = _pol_angle

#_mask = (_energy >= MIN_ENERGY)*(_energy <= MAX_ENERGY)
#_energy = _energy[_mask]
_pol_angle = _pol_angle[_mask]

print "pol angle", _pol_angle
print "energy",_energy

fmt = dict(xname='Energy', yname='Polarization angle [rad]')
pol_angle_spline = xInterpolatedUnivariateSpline(_energy, _pol_angle, k=1, **fmt)

#"""
#put together the flux for the source starting from the txt file
_energy, _flux = numpy.loadtxt(FLUX_FILE_PATH, unpack=True)
_mask = (_energy >= MIN_ENERGY)*(_energy <= MAX_ENERGY)
_energy = _energy[_mask]
_flux = _flux[_mask]
fmt = dict(xname='Energy', xunits='keV', yname='Flux',
               yunits='cm$^{-2}$ s$^{-1}$ keV$^{-1}$')

_energy_spectrum = xInterpolatedUnivariateSpline(_energy, _flux, **fmt)


ROI_MODEL = xROIModel(GRS1915_RA, GRS1915_DEC)


src = xPointSource('GRS1915_105', ROI_MODEL.ra, ROI_MODEL.dec, energy_spectrum,
                   polarization_degree, polarization_angle)

ROI_MODEL.add_source(src)



if __name__ == '__main__':
    print(ROI_MODEL)
    from ximpol.utils.matplotlib_ import pyplot as plt
    fig = plt.figure('Spectrum')
    _energy_spectrum.plot(show=False)
    fig = plt.figure('Polarization angle')
    pol_angle_spline.plot(show=False)
    fig = plt.figure('Polarization degree')
    pol_degree_spline.plot(show=False)
    plt.show()
